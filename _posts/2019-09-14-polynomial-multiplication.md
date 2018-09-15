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

School algorithm is very simple, but also very slow. It's time complexity is $O(n^2)$ .
There are much more efficient solutions and we want to explore them in this article.

After formalizing the problem, we'll see Karatsuba algorithm for multiplying polynomials in $O(n^{1.58})$ .
Next, we'll try to multiply polynomials using polynomial interpolation which will be the base for $O(n \log n)$ algorithm that utilizes Fast Fourirer Transform (FFT).

All these algorithms are [Divide and Conquer](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithms) algorithms, so this article can also serve as an introduction (through examples) with that problem solving approach.

## Problem definition
In the rest of the article we'll use following ways to note polynomial of degree $n$ (all are considered equal):

$p = p(x) = \sum\limits_{j=0}^n p_jx^j = p_0 + p_1x + p_2x^2 + .. + p_nx^n$

On the input we have two polynomials $a$ and $b$ of the same degree $n$, and the output is a polynomial $c$ of degree $2n$:

$c(x) = a(x)b(x) = \sum\limits_{i=0}^n \sum\limits_{j=0}^n a_ib_j x^{i+j}$

School algorithm written in C is simple and unambiguosly illustrates the problem:
{% highlight c %}

  void mnozi(int* a, int* b, int* c, int n) {
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
	* if we compute $z_0$ and $z_2$ directly, it turns out that we can compute $z_1$ as:
		* $z_1 = (p + q)(r + s) - pr - qs = (p + q)(r + s) - z_0 - z_2$

One multiplication is compensated by few additions (which are cheap) so now we need only **three multiplications** of half-sized polynomials.
We can multiply those by invoking the algorithm recursively and combine the results according to the last expression for $c$.

Time complexity of this algorithm is $O(n^{\log_2 3})$ which is close to, and usually written as, $O(n^{1.58})$.

If you are curious how to compute the complexity of this and similar algorithms in an elegant way, consult the [Master theorem (Wikipedia)](https://en.wikipedia.org/wiki/Master_theorem_(analysis_of_algorithms)).


## Multiplications using polynomial interpolation

Poznato je da $n+1$ parova $(x_0, y_0), .. , (x_n, y_n)$ takvih da $x_i \neq x_j$ za $i \neq j$ jednoznačno određuju polinom $p$ stupnja $n$ takav da $p(x_i) = y_i$ . Kažemo da iz evaluacija polinoma u $n+1$ točaka možemo rekonstruirati polinom.

Nama je u problemu množenja polinoma nepoznat polinom $c$ stupnja $2n$ . Ako pronađemo vrijednosti tog polinoma u $2n + 1$ točaka možemo iz njih rekonstruirati cijeli polinom!
Kako pronaći vrijednost polinoma $c$ u nekoj točki $x$ ? Budući da se radi o umnošku polinoma $a$ i $b$ dovoljno je evaluirati oba polinoma u točki $x$ i pomnožiti dobivene vrijednosti.

Algoritam se dakle sastoji od 3 koraka:

	* 1. Izaberi proizvoljno $2n+1$ različitih vrijednosti (točaka) $x_0, .. , x_{2n}$ i evaluiraj polinome $a$ i $b$ u njima.
	* 2. Pomnoži vrijednosti dobivene za svaku točku: $y_i = a(x_i) * b(x_i)$ .
	* 3. Interpolacijom izračunaj polinom $c$ iz skupa parova $(x_0, y_0), .. , (x_{2n}, y_{2n})$ .

Sve je to lijepo i krasno, no već prvi korak ne znamo napraviti u složenosti boljoj od $O(n^2)$ .
Netko je skužio da se evaluacija polinoma u hrpi točaka može napraviti brže ako su točke odabrane prema određenim pravilima. I to nije sve, treći korak se u tom slučaju isto može svesti na evaluaciju polinoma sličnog skupa točaka. Drugi korak je očito linearan, pa ovo zvuči obećavajuće.
Upoznajmo se najprije s brzim postupkom evaluacije polinoma u takvim točkama.

### Multipoint evaluation

Algoritam koji ćemo u ovom dijelu proučiti uzima polinom $p$ , prirodan broj $N$ i $w$ .
Algoritam vraća vrijednosti polinoma $p$ u $N$ točaka: $w^0$ , $w^1$ , .. , $w^{N-1}$ .
Zgodno je vrijednosti dobivene evaluacijom zapisati u obliku polinoma stupnja $N-1$ , tako da vrijednost evaluacije u točki $w^j$ stoji uz $x^j$ . U tom slučaju možemo reći da algoritam vraća novi polinom, i tada govorimo o transformaciji polinoma. **U ostatku ovog članka kad kažemo transformacija mislimo na ovaj algoritam.**

Matematički, transformacija koju radimo je:
$T_{N, w}(p) = \sum\limits_{j=0}^{N-1}p(w^j)x^j$

Dodat ćemo još jedan uvjet koji $w$ mora zadovoljavati (ne pitajte zašto, barem ne zasad): $w^N = 1$ . Ako ste upravo pomislili da uvjet nema previše smisla ne zaboravite da se nismo ograničili na realne brojeve kod odabira $w$ .

Kako brzo izračunati ovu transformaciju?
Upoznat ćemo se s varijantom **Cooley-Tukey** algoritma koja brzo izračunava ovu transformaciju kada je $N$ potencija broja $2$ .

Na ulazu su dakle $N=2^k$ , $w$ i polinom $p$ , a na izlazu je transformirani polinom $p' = T_{N, w}(p)$ .

Algoritam je rekurzivan:

	* temeljni slučaj je $N = 1$
		* $p' = p(w^0) = p(1) = p_0$
	* inače:
		* promotrimo neki član rezultata $p'_j$ :
			* $p'_j = p(w^j) = \sum\limits_{k=0}^{N-1}w^{jk}p_k$
		* rastavimo ga po parnosti potencija u $p$ , neka je $m = N/2$ :
			* $p'_j = \sum\limits_{k=0}^{m-1}w^{j2k}p_{2k} + \sum\limits_{k=0}^{m-1}w^{j(2k+1)}p_{2k+1}$
			* $p'_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}p_{2k} + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}p_{2k+1}$
		* uvedimo supstitucijske polinome $e$ i $o$ , koji sadrže koeficijente uz parne, odnosno neparne potencije iz $p$ :
			* $e_j = p_{2j}$
			* $o_j = p_{2j+1}$
		* napišimo sada $p'_j$ koristeći $o$ i $e$ :
			* $p'_j = \sum\limits_{k=0}^{m-1}(w^2)^{jk}e_k + w^j\sum\limits_{k=0}^{m-1}(w^2)^{jk}o_k$
			* $p'_j = e((w^2)^j) + w^jo((w^2)^j)$
		* ono što je sada potrebno prepoznati u posljednjem izrazu jesu izrazi za transformaciju polinoma $e$ i $o$ uz kvadrirani $w$
		* transformirajmo rekurzivno polinome $e$ i $o$ s kvadriranim $w$  (primjetimo $(w^2)^m = 1$ )
			* $e' = T_{m, w^2}(e)$
			* $o' = T_{m, w^2}(o)$
		* ako se sada vratimo na posljednji izraz za $p'$ očito je da za $0 \leq j \lt m$ vrijedi:
			* $p'_j = e'_j + w^jo'_j$
		* problem sa $j \ge m$ je što manje polinome $e$ i $o$ pri rekurzivnoj transformaciji nismo evaluirali u točki $(w^2)^j$ , no budući da je $(w^2)^m = 1$ možemo iskoristiti evaluaciju u točki $(w^2)^{j-m}$ tj. $e'_{j-m}$ odnosno $o'_{j-m}$ .

		* dakle:
			* $p'_j =
  				\begin{cases}
    				e'_j + w^jo'_j       & \quad \text{za } 0 \le j \lt m\\
    				e'_{j-m} + w^jo'_{j-m}  & \quad \text{za } m \le j \lt N \\
  					\end{cases}$


Složenost algoritma je $O(N\log N)$ .


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
