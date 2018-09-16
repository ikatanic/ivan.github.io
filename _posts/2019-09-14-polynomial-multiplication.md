---
layout: post
title:  "How to efficiently multiply polynomials"
date:   2018-09-13 08:58:09 +0100
categories: jekyll update
---

* TOC
{:toc}

## Motivation
Many problems can be reduced to a polynomial multiplication problem.
The most familiar example is the number multiplication.
Have you ever wondered why are you able to multiply 100000-digit numbers in Python in seconds?

School algorithm is very simple, but also very slow. Its time complexity is $O(n^2)$.
There are much more efficient algorithms and the goal of this article is to explore them.

After formalizing the problem, we'll see the Karatsuba algorithm for multiplying polynomials in $O(n^{1.58})$.
Next, we'll try to multiply polynomials using polynomial interpolation which will be the base for $O(n \log n)$ algorithms that utilize Fast Fourier Transform (FFT) and Number Theoretic Transform (NTT).

All these algorithms are [Divide and Conquer](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithms) algorithms, so this article can also serve as an introduction (through examples) with that problem-solving approach.

## Problem definition
In the rest of the article we'll use following ways to note polynomial of degree $n$ (all are considered equal):

$p = p(x) = \sum\limits_{j=0}^n p_jx^j = p_0 + p_1x + p_2x^2 + .. + p_nx^n$

At the input we have two polynomials $a$ and $b$ of same degree $n$, and the output is a polynomial $c$ of degree $2n$:

$c(x) = a(x)b(x) = \sum\limits_{i=0}^n \sum\limits_{j=0}^n a_ib_j x^{i+j}$

School algorithm written in C is simple and unambiguously illustrates the problem:
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
* base case is $n = 0$, then simply $c = a_0b_0$
* otherwise ($n > 0$):
	* split the input polynomials into halves of (almost) equal sizes:
		* $m = \lfloor n/2 \rfloor$
		* $a = px^m + q$
		* $b = rx^m + s$
	* write $c$ using new smaller polynomials $p, q, r, s$:
		* $c = ab = (px^m + q)(rx^m + s) = prx^{2m} + (ps + qr)x^m + qs$
	* introduce simple substitutions:
		* $z_2 = pr$
		* $z_1 = ps + qr$
		* $z_0 = qs$
	* now we have:
		* $c = z_2x^{2m} + z_1x^m + z_0$
	* we see that we can compute $z_0, z_1$ and $z_2$ with **four multiplications** of half-sized polynomials, but we can do better
	* if we compute $z_0$ and $z_2$ directly, it turns out that we can compute $z_1$ as follows:
		* $z_1 = (p + q)(r + s) - pr - qs = (p + q)(r + s) - z_0 - z_2$

One multiplication is compensated by few additions (which are cheap) so now we need only **three multiplications** of half-sized polynomials.
We can multiply those by invoking the algorithm recursively and combine the results according to the last expression for $c$.

Time complexity of this algorithm is $O(n^{\log_2 3})$ which is close to, and usually written as, $O(n^{1.58})$.

If you are curious about how to compute the complexity of this and similar algorithms in an elegant way, consult the [Master theorem (Wikipedia)](https://en.wikipedia.org/wiki/Master_theorem_(analysis_of_algorithms)).


## Multiplications using polynomial interpolation

It's known that $n+1$ pairs $(x_0, y_0), .. , (x_n, y_n)$ such that $x_i \neq x_j$ for $i \neq j$ uniquely determine a polynomial $p$ of degree $n$ such that $p(x_i) = y_i$. We say that from evaluations of a polynomial in $n+1$ points we can reconstruct (interpolate) the polynomial. Simplest intuition for this is that a polynomial has $n + 1$ unknown coefficients, and each evaluation is like an equation.

In the polynomial multiplication problem, we are looking for the polynomial $c$ of degree $2n$. What if we find values of $c$ in $2n + 1$ points? Then we can reconstruct the solution!

How to find the value of $c$ in a point $x$? As $c$ is the product of $a$ and $b$ it is enough to evaluate both polynomials in $x$ and multiply the results.

Now we are ready to sketch out the high level 3-step algorithm:

1. Choose some $2n+1$ different numbers (points) $x_0, .. , x_{2n}$ and evaluate polynomials $a$ and $b$ in them,
2. Multiply evaluations for each point: $y_i = a(x_i) b(x_i)$,
3. Interpolate $c$ from pairs $(x_0, y_0), .. , (x_{2n}, y_{2n})$.

This is great, but already the first step is difficult to do better than $O(n^2)$.
Lucky for us, someone smart noticed that evaluation of a polynomial in many points (multipoint evaluation) can be done more efficiently if the chosen set of points satisfies few rules.
And that's not everything, the third step, in that case, can also be reduced to the multipoint evaluation of a similar set of points.
The second step is clearly linear, so this sounds promising.

Let's first see how to do multipoint evaluation efficiently.

### Multipoint evaluation

The algorithm we are going to see takes a polynomial $p$, a natural number $N$ and a complex number $w$.
It returns evaluations of $p$ in $N$ points: $w^0, w^1,.. , w^{N-1}$.

It's convenient to write down these $N$ evaluations in the form of a degree $N-1$ polynomial, so that value of $p$ in $w^j$ stays next to $x^j$.
In that case, we can say that the algorithm returns a new polynomial, so we can also say the algorithm is a polynomial transformation.
It takes a polynomial, it returns a polynomial.
In the rest of the article, when we say **transformation** we mean this algorithm.

Formally, the transformation is:

$T_{N, w}(p) = \sum\limits_{j=0}^{N-1}p(w^j)x^j$

With an extra condition: $w^N = 1$ (later we'll see why).

How to compute this transformation efficiently?
We'll use a variant of **Cooley-Tukey** algorithm that computes this transformation when $N$ is a power of 2.


#### Cooley-Tukey algorithm for multipoint evaluation

Input: $N=2^k$, $w$ and polynomial $p$.

Output: transformed polynomial $p\prime = T_{N, w}(p)$ .

Again, we solve recursively:
* base case is $N = 1$
	* $p\prime = p(w^0) = p(1) = p_0$
* otherwise $(N > 1)$:
	* let's look at one element of the result $p\prime_j$:
		* $p\prime_j = p(w^j) = \sum\limits_{k=0}^{N-1}w^{jk}p_k$
	* split on parity of powers in $p$, and let $m = N/2$:
		* $p\prime_j = \sum\limits_{k=0}^{m-1}w^{j2k}p_{2k} + \sum\limits_{k=0}^{m-1}w^{j(2k+1)}p_{2k+1}$
		* $p\prime_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}p_{2k} + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}p_{2k+1}$
	* introduce substitutions $e$ and $o$, for coefficients next to even and odd powers in $p$, respectively:
		* $e_j = p_{2j}$
		* $o_j = p_{2j+1}$
	* express $p\prime_j$ using $e$ i $o$ :
		* $p\prime_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}e_k + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}o_k$
		* $p\prime_j = e((w^2)^j) + w^jo((w^2)^j)$
	* here we have to recognize transformation expressions of polynomials $e$ and $o$ with squared $w$
	* transform recursively $e$ and $o$ with $w^2$ (note that $(w^2)^m = 1$ so we are allowed to do it):
		* $e\prime = T_{m, w^2}(e)$
		* $o\prime = T_{m, w^2}(o)$
	* now it's clear that for $0 \leq j \lt m$:
		* $p\prime_j = e\prime_j + w^jo\prime_j$
	* what about $j \ge m$? it might seem impossible because we haven't evaluated $e$ and $o$ in points $(w^2)^j$,
	but since $(w^2)^m = 1$ we can use evaluations in $(w^2)^{j-m}$, that is, $e\prime_{j-m}$ or $o\prime_{j-m}$.
	* to summarize:
		* $$p\prime_j =
 				\begin{cases}
   				e\prime_j + w^jo\prime_j          & \quad \text{for } 0 \le j \lt m \\
   				e\prime_{j-m} + w^jo\prime_{j-m}  & \quad \text{for } m \le j \lt N \\
				\end{cases}$$

The complexity of this multipoint evaluation algorithm is $O(N\log N)$.

Cool, we can evaluate efficiently. But what about the interpolation?
We can look at the interpolation as an inverse of our transformation.
It turns out it is the same as the transformation with inverted $w$ and divided by $N$:

$$T^{-1}_{N, w}(p\prime) = \frac{1}{N} T_{N, w^{-1}}(p\prime) = \frac{1}{N} \sum\limits_{j=0}^{N-1}p\prime(w^{-j})x^j$$

The proof will be left as an exercise to the reader :).

So basically we can use the same multipoint evaluation algorithm to do the interpolation too!

The only thing we are missing is a value of $w$...

### Multiplication using Fast Fourier Transform (FFT)

Let's use $w = e^{i2\pi /N}$.

In that case, our transform is also known as the Discrete Fourier Transform (DFT).

If we now go back to the original multiplication problem, we can summarize the complete algorithm as follows:

$$c = ab = DFT^{-1}(DFT(a).DFT(b))$$

where:
* $N = $ smallest power of two greater than $2n$,
* $DFT(p) = T_{N, e^{i2\pi /N}}(p)$,
* $DFT^{-1}(p) = \frac{1}{N} T_{N, e^{-i2\pi /N}}(p)$,
* $.$ - element-wise multiplication of polynomials, that is: $(r.s)_i = r_is_i$.

This and all other fast implementation of DFT are called Fast Fourier Transformations (FFT)	.

The complexity of this polynomial multiplication algorithm is  $O(n \log n)$.

### Multiplication using Number Theoretic Transform (NTT)

A disadvantage of DFT in the context of implementation can be the fact that it uses real numbers.
If we work with polynomials over finite fields we may go around it using Number Theoretic Transform (NTT).
Let's say we are working with a finite field $GP(p)$, that is, all the operations on numbers are done modulo $p$, where $p$ is prime.

Similarly to DFT, the product of polynomials $a$ and $b$ of degree $n$ can we written as:

$$c = ab = NTT^{-1}(NTT(a).NTT(b))$$

where:
* $N = $ smallest divisor of $p-1$ that is greater than $2n$,
* $g = $ a primitive root modulo $p$  [What is primitive root and how to find one?](https://en.wikipedia.org/wiki/Primitive_root_modulo_n)
* $NTT(p) = T_{N, g^{(p-1)/N}}(p)$,
* $NTT^{-1}(p) = \frac{1}{N} T_{N, g^{-(p-1)/N}}(p)$,
* $.$ - element-wise multiplication of polynomials, that is: $(r.s)_i = r_is_i$.

Note that $g^{(p-1)/N}$ has similar properties as our old friend $w$. So we can simply use the
same transformation algorithm as before, the main difference being that all the computations will be done using integers and modulo $p$.

There is a catch though. $N$ must divide $p-1$, and with the old requirement of $N$ being a power of two, finding an appropriate $N$ can become impossible.
To summarize, to apply this algorithm, $p$ must be a prime of form $k2^l + 1$.

However, it is sometimes possible to modify the Cooley-Tukey to support other values of $p$, by not always splitting the problem into two equal parts.
Details are out of the scope of this article, but feel free to experiment.


## Conclusion
We started with the polynomial multiplication problem but we also learned how to do FFT efficiently. FFT, on the other hand, is used everywhere (signal processing for example).

FFT/NTT-based multiplications are definitely a better choice than Karatsuba when talking about speed, but FFT can have
precision problems while NTT might not be applicable.
Another random benefit of Karatsuba algorithm is that it can be used even if the division is not defined (modulo 24, for example).

## Few practical tips
1. NTT-based multiplication can be used even if we don't want the result modulo prime $p$. It is sufficient to use any $p$ larger than any value in our result.
If we are constrained by sizes of integer data types, we can do the multiplication modulo few different primes. They can be small, but their product must be larger
than any value in our result. The results can then be combined into non-modulo result using [Chinese Remained Theorem](https://en.wikipedia.org/wiki/Chinese_remainder_theorem).

2. When implementing Cooley-Tukey to do FFT be aware of the form of $w$ and compute its powers using [De Moivre's formula](https://en.wikipedia.org/wiki/De_Moivre%27s_formula).

3. To combat FFTs precision issues with big numbers we can split the polynomial's coefficients into few smaller values. For example
each coefficient $a_i$ can be expressed as $a_i = x_iK + y_i$. Where $K$ is some constant, for example $K = \sqrt{max(a_i)}$.
If we do a similar split for the other polynomial, we can express the product as a weighted sum of 4 smaller (in terms of values) products.

## Test your understanding
This article was originally written as a competitive programming tutorial. So naturally, I've collected a few problems available online, ranging from ones where you just have to implement one of the algorithms, to ones where you have to modify the algorithms to adapt them to a new situation:

* [SPOJ VFMUL](http://www.spoj.com/problems/VFMUL/)
* [SPOJ LPRIME](http://www.spoj.com/problems/LPRIME/)
* [SPOJ MAXMATCH](http://www.spoj.com/problems/MAXMATCH/)
* [AGC 5 - problem F](https://agc005.contest.atcoder.jp/tasks/agc005_f)
* [ICPC World Finals 2015 - problem J](http://icpc.baylor.edu/download/worldfinals/problems/icpc2015.pdf)
* [ASC 46 - problem D](http://codeforces.com/gym/100524)
