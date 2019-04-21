---
layout: post
title:  "Fast polynomial multiplication"
date:   2018-09-22 08:58:09 +0100
categories: cp math algorithms
---

* TOC
{:toc}

## Motivation
Many computer science problems can be reduced to a polynomial multiplication.
The most familiar example is the integer multiplication.
Have you ever wondered why multiplying 100000-digit numbers in Python takes seconds and not minutes?

The school algorithm is very simple, but also very slow. Its time complexity is $O(n^2)$.
There are much more efficient algorithms and the goal of this article is to explore some of them.

After formalizing the problem, we'll have a look the Karatsuba algorithm for multiplying polynomials in $O(n^{1.58})$ time.
Next, we'll try to multiply polynomials using polynomial interpolation which will be the basis of $O(n \log n)$ algorithms that utilize Fast Fourier Transform (FFT) and Number Theoretic Transform (NTT).

All these algorithms use the [Divide and Conquer](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithms) paradigm, so this article can also serve as an advanced introduction (through examples) to that problem-solving approach.

## Problem definition
In the rest of the article we'll denote a polynomial of degree $n$ using one of the following:

$$p = p(x) = \sum\limits_{j=0}^n p_jx^j = p_0 + p_1x + p_2x^2 + .. + p_nx^n$$

The input consists of two polynomials $a$ and $b$ of the same degree $n$ (if that's not the case we can always pad with zeroes), and the output is a polynomial $c$ of degree $2n$:

$$c(x) = a(x)b(x) = \sum\limits_{i=0}^n \sum\limits_{j=0}^n a_ib_j x^{i+j}$$

The school algorithm written in C is simple and unambiguously illustrates the problem:
{% highlight c %}

  void mul(int* a, int* b, int* c, int n) {
    for (int i = 0; i <= 2*n; ++i) {
      c[i] = 0;
    }
    for (int i = 0; i <= n; ++i) {
      for (int j = 0; j <= n; ++j) {
        c[i+j] += a[i] * b[j];
      }
    }
  }
{% endhighlight %}




## Karatsuba algorithm
We will solve the problem recursively:
* the base case is $n = 0$, then simply $c = a_0b_0$
* otherwise ($n > 0$):
	* split each of the input polynomials into two smaller polynomials:
		* $m = \lceil n/2 \rceil$
		* $a = px^m + q$
		* $b = rx^m + s$
		* for example, if $a = x^3 - 2x^2 + x + 4$, then $m=2$, $q=x + 4$ and $p=x - 2$
	* write $c$ using new polynomials $p, q, r, s$:
		* $c = ab = (px^m + q)(rx^m + s) = prx^{2m} + (ps + qr)x^m + qs$
	* introduce simple substitutions:
		* $z_2 = pr$
		* $z_1 = ps + qr$
		* $z_0 = qs$
	* now we have:
		* $c = z_2x^{2m} + z_1x^m + z_0$
	* we see that we can compute $z_0, z_1$ and $z_2$ with **four multiplications** of half-sized polynomials, but we can do better
	* if we compute $z_0$ and $z_2$ first, it turns out we can compute $z_1$ as follows:
		* $z_1 = (p + q)(r + s) - pr - qs = (p + q)(r + s) - z_0 - z_2$

One multiplication is substituted by several additions (which are cheap) so now we need only **three multiplications** of half-degreed polynomials.
We can multiply those by invoking the algorithm recursively.

Time complexity of this algorithm is $O(n^{\log_2 3})$ which is close to, and is usually written as, $O(n^{1.58})$.

If you are curious about how to compute the complexity of this and similar algorithms in an elegant way, see the [Master theorem (Wikipedia)](https://en.wikipedia.org/wiki/Master_theorem_(analysis_of_algorithms)).


## Multiplications using polynomial interpolation

It is known that $n+1$ pairs $(x_0, y_0), .. , (x_n, y_n)$ where $x_i \neq x_j$ for $i \neq j$ uniquely determine a polynomial $p$ of degree $n$ where $p(x_i) = y_i$. We say that from evaluations of a polynomial in $n+1$ points we can reconstruct (interpolate) the polynomial. 
Simplest intuition for this is that the polynomial has $n + 1$ unknown coefficients, and each evaluation is an equation.

In the polynomial multiplication problem, we are looking for the polynomial $c$ of degree $2n$. What if we find values of $c$ in $2n + 1$ points? Then we should be able to reconstruct the solution!

How to find the value of $c$ in a point $x$? Since $c(x) = a(x)b(x)$, it suffices to evaluate both polynomials in $x$ and multiply the results.

Now we are ready to sketch out the high-level, three-step algorithm:

1. Choose $2n+1$ different numbers (points) $x_0, .. , x_{2n}$ and evaluate polynomials $a$ and $b$ in them,
2. Multiply evaluations for each point: $y_i = a(x_i) b(x_i)$,
3. Interpolate $c$ from pairs $(x_0, y_0), .. , (x_{2n}, y_{2n})$.

This is great, but already the first step is difficult to perform in better than $O(n^2)$ time.
Lucky for us, someone noticed that evaluation of a polynomial in many points (multipoint evaluation) can be done more efficiently if the chosen set of points satisfies several rules.
Furthermore, the third step can then also be reduced to the multipoint evaluation of a similar set of points.
The second step is clearly linear, so this sounds promising.

Let's first see how to do the multipoint evaluation efficiently.

### Multipoint evaluation

The algorithm we are going to see takes a polynomial $p$, a natural number $N$ and a complex number $w$.
It returns evaluations of $p$ in $N$ points: $w^0, w^1,.. , w^{N-1}$.

It's convenient to write down these $N$ evaluations in the form of a degree $N-1$ polynomial which is their
[generating function](https://en.wikipedia.org/wiki/Generating_function).
In that case, we can say that the algorithm returns a new polynomial, so we can also say the algorithm is a polynomial transformation.
It takes a polynomial and returns a polynomial.
In the rest of the article, when we say **transformation** we mean this transform.

Formally, the transformation is:

$$T_{N, w}(p) = \sum\limits_{j=0}^{N-1}p(w^j)x^j$$

In general, it is hard to compute this transformation efficiently. But not if $w$ is chosen carefully.
As we'll see in the following section, the transformation can be computed efficiently
if $w$ is an $N$th [root of unity](https://en.wikipedia.org/wiki/Root_of_unity), that is, if $w^N=1$.

We'll use a variant of **Cooley-Tukey** algorithm that computes this transformation when $N$ is a power of 2.


#### Cooley-Tukey algorithm for multipoint evaluation

Input: $N=2^q$, $w$ such that $w^N=1$ and a polynomial $p$.

Output: transformed polynomial $p\prime = T_{N, w}(p)$ .

Again, we solve recursively:
* the base case is $N = 1$
	* $p\prime = p(w^0) = p(1) = p_0$
* otherwise $(N > 1)$:
	* let's look at one element of the result $p\prime_j$:
		* $p\prime_j = p(w^j) = \sum\limits_{k=0}^{N-1}w^{jk}p_k$
	* split on parity of powers in $p$, and let $m = \frac{N}{2}$:
		* $p\prime_j = \sum\limits_{k=0}^{m-1}w^{j2k}p_{2k} + \sum\limits_{k=0}^{m-1}w^{j(2k+1)}p_{2k+1}$
		* $p\prime_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}p_{2k} + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}p_{2k+1}$
	* introduce substitutions $e$ and $o$, for coefficients next to even and odd powers in $p$, respectively:
		* $e = \sum\limits_{j=0}^{m-1}p_{2j}x^j$
		* $o = \sum\limits_{j=0}^{m-1}p_{2j+1}x^j$
	* express $p\prime_j$ using $e$ i $o$ :
		* $p\prime_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}e_k + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}o_k$
		* $p\prime_j = e((w^2)^j) + w^jo((w^2)^j)$
	* here we have to recognize transformation expressions of polynomials $e$ and $o$ with squared $w$
	* transform recursively $e$ and $o$ with $w^2$ (note that $(w^2)^m = 1$):
		* $e\prime = T_{m, w^2}(e)$
		* $o\prime = T_{m, w^2}(o)$
	* now it's clear that for $0 \leq j \lt m$:
		* $p\prime_j = e\prime_j + w^jo\prime_j$
	* what about $j \ge m$? It might seem impossible because we haven't evaluated $e$ and $o$ in points $(w^2)^j$,
	but since $(w^2)^m = 1$ we can use the evaluations in $(w^2)^{j-m} = (w^2)^j$.
	* to summarize:
		* $$p\prime_j =
 				\begin{cases}
   				e\prime_j + w^jo\prime_j          & \quad \text{for } 0 \le j \lt m \\
   				e\prime_{j-m} + w^jo\prime_{j-m}  & \quad \text{for } m \le j \lt N \\
				\end{cases}$$

Implementing it recursively will give us an algorithm of $O(N\log N)$ time complexity.

Cool, we can evaluate efficiently. But what about the interpolation?
We can look at the interpolation as an inverse of our transformation.
It turns out it is the same as the transformation with inverted $w$ and divided by $N$:

$$T^{-1}_{N, w}(p\prime) = \frac{1}{N} T_{N, w^{-1}}(p\prime) = \frac{1}{N} \sum\limits_{j=0}^{N-1}p\prime(w^{-j})x^j$$

The proof is left as an exercise to the reader :).

What it means is that we can use the same multipoint evaluation algorithm to do the interpolation too!

The only thing missing is the value of $w$...

### Multiplication using Fast Fourier Transform (FFT)

Let's use $w = e^{i2\pi /N}$.

In that case, our transform is also known as the Discrete Fourier Transform (DFT).

If we now go back to the original multiplication problem, we can summarize the complete algorithm as follows:

$$c = ab = DFT^{-1}(DFT(a) \cdot DFT(b))$$

where:
* $DFT(p) = T_{N, e^{i2\pi /N}}(p)$,
* $DFT^{-1}(p) = \frac{1}{N} T_{N, e^{-i2\pi /N}}(p)$,
* $N = $ smallest power of two greater than $2n$ and
* $\cdot$ is element-wise multiplication operator on polynomials, i.e.: $(r \cdot s)_i = r_is_i$.

We can use the Cooley-Tukey algorithm to compute the DFT efficiently.
Such and all other fast implementations of DFT are called Fast Fourier Transformations (FFT).

The complexity of this polynomial multiplication algorithm is  $O(n \log n)$.

### Multiplication using Number Theoretic Transform (NTT)

A disadvantage of DFT in the context of implementation can be the fact that it uses complex numbers.
If we work with polynomials over finite fields we may go around it using Number Theoretic Transform (NTT).
Let's say we are working with a finite field $GP(p)$, that is, all the operations on numbers are done modulo $p$, where $p$ is prime.

Similarly to DFT, the product of polynomials $a$ and $b$ of degree $n$ can we written as:

$$c = ab = NTT^{-1}(NTT(a) \cdot NTT(b))$$

where:
* $NTT(p) = T_{N, g^{(p-1)/N}}(p)$,
* $NTT^{-1}(p) = \frac{1}{N} T_{N, g^{-(p-1)/N}}(p)$,
* $g = $ a primitive root modulo $p$  ([What is a primitive root and how to find one?](https://en.wikipedia.org/wiki/Primitive_root_modulo_n)),
* $N = $ smallest divisor of $p-1$ that is greater than $2n$ and
* $\cdot$ is element-wise multiplication operator on polynomials, i.e.: $(r \cdot s)_i = r_is_i$.

Note that $(g^{(p-1)/N)})^N=1$. Hence we can use the same transformation algorithm as before,
the main difference being that all the computations will be done using integers and modulo $p$.

There is a catch though. $N$ must divide $p-1$, and with the old requirement of $N$ being a power of two, finding an appropriate $N$ can be an impossible task.

To summarize, to apply this algorithm, $p$ must be a prime of form $k2^q + 1$ and $N$ must be a power of two and at most $2^q$.

However, it is sometimes possible to modify the Cooley-Tukey to support other values of $p$, by not always splitting the problem into two equal parts.
Details are out of the scope of this article, but feel free to experiment.


## Conclusion
We started with the polynomial multiplication problem but we also learned how to do FFT efficiently. FFT, on the other hand, is used everywhere (for example, processing of various kinds of signals). Some big-integer libraries still use the Karatsuba algorithm, while others have opted for FFT or even fancier algorithms.

FFT/NTT-based multiplications are definitely a better choice than Karatsuba when talking about speed, but FFT can have precision problems while NTT might not be applicable.
Another random benefit of the Karatsuba algorithm is that it can be used even when the division is not defined (modulo 24, for example).

## Few practical tips
1. NTT-based multiplication can be used even if we don't want the result modulo prime $p$. It is sufficient to use any $p$ larger than any value in our result.
If we are constrained by sizes of integer data types, we can do the multiplication modulo few different primes. They can be small, but their product must be larger
than any value in our result. The results can then be combined into non-modulo result using the [Chinese Remained Theorem](https://en.wikipedia.org/wiki/Chinese_remainder_theorem).

2. When implementing Cooley-Tukey to do FFT be aware of the form of $w$ and compute its powers using [De Moivre's formula](https://en.wikipedia.org/wiki/De_Moivre%27s_formula).

3. To combat FFTs precision issues with big numbers we can split the polynomial's coefficients into few smaller values. For example
each coefficient $a_i$ could be expressed as $a_i = x_iK + y_i$. Where $K$ is some constant, for example $K = \sqrt{max(a_i)}$.
If we do a similar split for the other polynomial, we can express the product as a weighted sum of 4 smaller (in terms of values) products.

4. Karatsuba algorithm will be faster than more sophisticated algorithms on smaller degree polynomials as it has a relatively low overhead factor. Similarly, the school algorithm, having the lowest overhead factor, will be the fastest for very small polynomials. So when the performance is critical, multiple approaches should be combined.

## Test your understanding
This article was originally written as a competitive programming tutorial. So naturally, I've collected a few interesting problems available online, ranging from the ones where you just have to implement one of these algorithms, to the ones where you have to modify the algorithms and adapt them:

* [SPOJ VFMUL](http://www.spoj.com/problems/VFMUL/)
* [SPOJ LPRIME](http://www.spoj.com/problems/LPRIME/)
* [SPOJ MAXMATCH](http://www.spoj.com/problems/MAXMATCH/)
* [AGC 5 - problem F](https://agc005.contest.atcoder.jp/tasks/agc005_f)
* [ICPC World Finals 2015 - problem J](http://icpc.baylor.edu/download/worldfinals/problems/icpc2015.pdf)
* [ASC 46 - problem D](http://codeforces.com/gym/100524)
