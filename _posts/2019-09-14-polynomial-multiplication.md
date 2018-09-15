---
layout: post
title:  "How to efficiently multiply polynomials"
date:   2018-09-13 08:58:09 +0100
categories: jekyll update
---

* TOC
{:toc}

## Motivation
Many problems can be reduced to polynomial multiplication problem.
Most familiar example is the number multiplication.
Have you ever wondered why are you able to multiply million-digit numbers in Python in seconds?

School algorithm is very simple, but also very slow. It's time complexity is $O(n^2)$.
There are much more efficient algorithms and the goal of this article is to explore them.

After formalizing the problem, we'll see Karatsuba algorithm for multiplying polynomials in $O(n^{1.58})$ .
Next, we'll try to multiply polynomials using polynomial interpolation which will be the base for $O(n \log n)$ algorithm that utilizes Fast Fourirer Transform (FFT).

All these algorithms are [Divide and Conquer](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithms) algorithms, so this article can also serve as an introduction (through examples) with that problem solving approach.

## Problem definition
In the rest of the article we'll use following ways to note polynomial of degree $n$ (all are considered equal):

$p = p(x) = \sum\limits_{j=0}^n p_jx^j = p_0 + p_1x + p_2x^2 + .. + p_nx^n$

At the input we have two polynomials $a$ and $b$ of same degree $n$, and the output is a polynomial $c$ of degree $2n$:

$c(x) = a(x)b(x) = \sum\limits_{i=0}^n \sum\limits_{j=0}^n a_ib_j x^{i+j}$

School algorithm written in C is simple and unambiguosly illustrates the problem:
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

If you are curious how to compute the complexity of this and similar algorithms in an elegant way, consult the [Master theorem (Wikipedia)](https://en.wikipedia.org/wiki/Master_theorem_(analysis_of_algorithms)).


## Multiplications using polynomial interpolation

It's known that $n+1$ pairs $(x_0, y_0), .. , (x_n, y_n)$ such that $x_i \neq x_j$ for $i \neq j$ uniquely determine a polynomial $p$ of degree $n$ such that $p(x_i) = y_i$. We say that from evaluations of a polynomial in $n+1$ points we can reconstruct (interpolate) the polynomial. Simplest intuition for this is that a polynomial has $n + 1$ unknown coefficients, and each evaluation is like an equation.

In the polynomial multiplication problem we are looking for the polynomial $c$ of degree $2n$. What if we find values of $c$ in $2n + 1$ points? Then we can reconstruct the whole solution!

How to find value of $c$ in a point $x$? As $c$ is the product of $a$ and $b$ it is enough to evaluate both polynomials in $x$ and multiply the results.

Now we are ready to sketch out the high level algorithm in 3 steps:

1. Choose arbitrarily $2n+1$ different numbers (points) $x_0, .. , x_{2n}$ and evaluate polynomials $a$ and $b$ in them.
2. Multiply evaluations for each point: $y_i = a(x_i) * b(x_i)$ .
3. Interpolate $c$ from pairs $(x_0, y_0), .. , (x_{2n}, y_{2n})$ .

This is great, but already first step is difficult to do better than $O(n^2)$.
Lucky for us, someone smart noticed that evaluation of a polynomial in many points (multipoint evaluation) can be done more efficiently if the chosen set of points satisfies some rules.
And that's not everything, third step in that case can also be reduced to multipoint evaluation of a similar set of points.
Second step is clearly linear, so this sounds promising.

Let's see first how to do multipoint evaluation efficiently.

### Multipoint evaluation

Algorithm we are going to see takes a polynomial $p$, natural number $N$ and a complex number $w$.
It returns evaluations of $p$ in $N$ points: $w^0$ , $w^1$ , .. , $w^{N-1}$.

It's convenient to write down these $N$ evaluations in the form of a degree $N-1$ polynomial, so that value of $p$ in $w^j$ stays next to $x^j$.
In that case, we can say that the algorithm returns a new polynomial, so we can also say the algorithm is a transformation of a polynomial.
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
	* write $p\prime_j$ using $e$ i $o$ :
		* $p\prime_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}e_k + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}o_k$
		* $p\prime_j = e((w^2)^j) + w^jo((w^2)^j)$
	* here we have to recognize transformation expressions of polynomials $e$ and $o$ with squared $w$
	* transform recursively $e$ and $o$ with $w^2$ (note that $(w^2)^m = 1$):
		* $e\prime = T_{m, w^2}(e)$
		* $o\prime = T_{m, w^2}(o)$
	* now it's clear that for $0 \leq j \lt m$:
		* $p\prime_j = e\prime_j + w^jo\prime_j$
	* what about $j \ge m$? it might seem impossible because we haven't evaluated $e$ and $o$ in points $(w^2)^j$,
	but since $(w^2)^m = 1$ we can use evaluations in $(w^2)^{j-m}$, that is, $e\prime_{j-m}$ or $o\prime_{j-m}$.
	* to summarize:
		* $p\prime_j =
 				\begin{cases}
   				e\prime_j + w^jo\prime_j          & \quad \text{for } 0 \le j \lt m \\
   				e\prime_{j-m} + w^jo\prime_{j-m}  & \quad \text{for } m \le j \lt N \\
				\end{cases}$

Complexity of the algorithm is $O(N\log N)$ .


### Multiplication using Fast Fourier Transform (FFT)

Ovaj algoritam nam kaže da umnožak polinoma $a$ i $b$ možemo zapisati kao:
$c = ab = DFT^{-1}(DFT(a).DFT(b))$

gdje je :

	* $N = 2n + 1$ ,
	* $DFT(p) = T_{N, e^{i2\pi /N}}(p)$ ,
	* $DFT^{-1}(p) = \frac{1}{N} T_{N, e^{-i2\pi /N}}(p)$ ,
	* $.$ - binarna operacija nad polinomima, "množenje i-tog člana s i-tim": $(r.s)_i = r_is_i$ .


Primjetite sličnosti s već spomenutim množenjem polinoma uporabom interpolacije.


Primjetite i da su vrijednosti $w$ ovdje kompleksni brojevi te da vrijedi $w^N = 1$ .
Možemo primjeniti algoritam iz prethodnog poglavlja za brzo izračunavanje $DFT$ transformacije.
Činjenica da transformaciju znamo brzo izračunati samo za slučajeve kad je $N$ potencija broja 2 ne predstavlja veliki problem, dovoljno je za $N$ uzeti prvu sljedeću potenciju broja 2 koja nije manja od $2n+1$ .
Sve ovakve i slične brze implementacije DFT-a nazivaju se zajedničkim imenom brze Fourieove transformacije (FFT).

Valja napomenuti da pri implementaciji ovog algoritma i transformacije treba biti svjestan oblika broja $w$ i ne računati direktno njegove potencije (kako bi se očuvala preciznost) već koristiti Eulerov identitet koji kaže da za bilo koji realan $x$ vrijedi:
$e^{ix} = \cos x + i \sin x$ .

Složenost ovog algoritma množenja polinoma je $O(n \log n)$ .

### Multiplication using Number Theoretic Transform (NTT)

Nedostatak DFT-a je (barem za nas) to što koristi realne odnosno kompleksne brojeve.
Ako radimo s polinomima nad konačnim poljima to ponekad možemo izbjeći. Pretpostavimo da radimo s konačnim poljem $GF(p)$ tj. sve operacije nad brojevima provodimo modulo $p$ gdje je $p$ prost broj.

Analogno DFT-u, umnožak polinoma $a$ i $b$ stupnja $n$ možemo pisati kao:
$c = ab = NTT^{-1}(NTT(a).NTT(b))$

gdje je :

	* $N$ - najmanji djelitelj od $p-1$ koji nije manji od $2n+1$ ,
	* $g$ - bilo koji primitivni korijen modulo $p$  ([[https://en.wikipedia.org/wiki/Primitive_root_modulo_n|Što je primitivni korijen i kako ga pronaći?]]),
	* $NTT(p) = T_{N, g^{(p-1)/N}}(p)$ ,
	* $NTT^{-1}(p) = \frac{1}{N} T_{N, g^{-(p-1)/N}}(p)$ ,
	* $.$ - binarna operacija nad polinomima, "množenje i-tog člana s i-tim": $(r.s)_i = r_is_i$ .

Za izračunavanje transformacije možemo iskoristiti već opisani algoritam, s tim da sve operacije provodimo modulo $p$ . Primjetite da smo sada dobili dodatan uvjet za $N$ , on mora biti djelitelj od $p-1$ . Uz stare zahtjeve da $N$ nije manji od $2n+1$ te da mora biti potencija broja 2, pronalaženje pravog $N$ može postati nemoguće. To je upravo i glavna mana NTT-a, jer je sada jedini slučaj gdje ga možemo iskoristiti slučaj kada je $p$ prost broj oblika $k2^n + 1$ .
Zapravo, moguće je i češće, ali to zahtjeva modificiranje opisanog algoritma za izračunavanje transformacije, glavna ideja ostaje ista ali ulazni niz ne rastavljamo konstantno na dva jednaka dijela. Oni koje to više zanima neka slobodno eksperimentiraju.


## Conclusion

Pošli smo od problema množenja polinoma, no treba napomenuti da se FFT koristi ponajviše u računalnom procesiranju raznih signala (slike, zvuka..).
Kad govorimo o množenju polinoma, klasična primjena je množenje velikih brojeva (zapisujemo ih kao polinome u bazi 10 ili nekoj većoj/manjoj bazi).

FFT je svakako bolji izbor od Karatsubinog algoritma kad govorimo o brzini, no nedostatak mu je preciznost jer smo primorani koristiti realni tip podatka pri implementaciji.
Uz preciznost, prednost Karatsubinog algoritma je i mogućnost množenja polinoma čak i kad dijeljenje nije definirano (npr. modulo 24).


Trebamo napomenuti da se NTT može iskoristiti čak i kad ne želimo eksplicitno rezultat modulo prost $p$ . Dovoljno je uzeti $p$ veći od najveće vrijednosti članova rezultantnog polinoma i pokrenuti algoritam. Nekad to nije praktično jer smo obično ograničeni 32-bitnim ili 64-bitnim cjelobrojnim tipom podatka. Tada možemo nekoliko puta izvršiti množenje, svaki puta modulo različit prost broj $p$ . Odabrani prosti brojevi mogu biti proizvoljno mali ali njihov umnožak $u$ mora biti veći od najveće vrijednosti članova rezultatntnog polinoma. Vrijednost člana rezultata modulo $u$  (tj. pravu vrijednost rezultata) možemo zatim izračunati koristeći Kineski teorem o ostacima (CRT).

## Test your understanding

http://www.spoj.com/problems/VFMUL/
http://www.spoj.com/problems/LPRIME/
http://www.spoj.com/problems/MAXMATCH/
http://codeforces.com/gym/100524 (problem D)
http://icpc.baylor.edu/download/worldfinals/problems/icpc2015.pdf (problem J)
