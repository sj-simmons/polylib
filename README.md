# crypto senior sem
### From cryptography to cryptocurrencies

This course is about the state-of-the-art in modern cryptography &mdash; and,
therefore, the primary mathematical apparatus that has, of late, been artfully
fashioned into cryptocurrencies.  Cryptography is of course pervasive and
ubiquitous in the technological age, so that this course is ultimately about
certain cutting-edge information-theoretic advancements that are reshaping our
world.  In addition to crytographic primitives, their prevailing applications,
and underlying mathematics, topics include applications of zero-knowledge proofs
and their role in present-day FinTech.

Important: as far as speculating in cryptocurrencies, we think that you are likely
better off not doing so. Whether or not you take that advice, you might
do well to look under the hood of the new-fangled tech that is built atop a rather
striking assemblage of pure and applied Math and CS.

TLDR; this is not Sh*tcoin 101 (which can be found on Youtube, TikTok and
elsewhere); rather, this course is about new takes on exceptionally powerful
mathematics. And remember always to please use your math skills for good,
not for ill.

## Topics

### Getting started

* First install [polylib](https://github.com/sj-simmons/polylib).
* Then install [numlib](https://github.com/sj-simmons/numlib).

### Cryptography

#### Large primes

We will need some prime numbers that are fairly large  &mdash; 200 bits, say, for now; so primes of size
around <img alt="$2^{200}.$" src="svgs/8ead495e1e0f713f33e7a75e00656907.svg" valign=0.0px width="33.264976799999985pt" height="13.380876299999999pt"/>

One way to generate such a prime would be to iterate through numbers larger that <img alt="$2^{200}$" src="svgs/acc9da0b19643e432619eb386d21261f.svg" valign=0.0px width="27.87685064999999pt" height="13.380876299999999pt"/> until we
find one that is prime.  Alternatively, one could randomly generate sequences of zeros and ones of
length 200 and check if the corresponding decimal number is prime.  Python has a built-in function that
generates an integer from random bits:

```python
import random
decimal = random.getrandbits(200)
```
Of course, depending on whether the most significant random bit was zero or one, we might get a number somewhat
less than <img alt="$2^{200};$" src="svgs/e46661a0be2003a6fba87fa4a557cf10.svg" valign=-3.1963503000000086px width="33.264976799999985pt" height="16.5772266pt"/> so let us set the most significant bit to one and, while we are at it, set
the least significant bit also equal to 1 since primes beyond 2 must be odd:
```python
decimal |= (1 << numbits - 1) | 1
```
The variable **decimal** is now an integer whose binary representation has length 200 and both begins
and ends with 1; i.e., **decimal** is a random (depending on the robustness of **getrandbits**)
integer strictly between <img alt="$2^{200}$" src="svgs/acc9da0b19643e432619eb386d21261f.svg" valign=0.0px width="27.87685064999999pt" height="13.380876299999999pt"/> and <img alt="$2^{201}.$" src="svgs/26851296087784570cb81e1034cdb3f0.svg" valign=0.0px width="33.264976799999985pt" height="13.380876299999999pt"/>

Beyond using the fact that a prime larger than 2 must be odd, there are various other quick ways
to test whether a candidate odd integer <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> is *likely* a prime.  These include
[Fermat's primality test](https://en.wikipedia.org/wiki/Fermat_primality_test) which checks to see
if <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> *acts like* a prime: namely, whether <img alt="$a^{n-1} \equiv 1 \mod n$" src="svgs/d366fecb67ec763d517425469c45fd90.svg" valign=0.0px width="122.41226084999998pt" height="13.380876299999999pt"/> for <img alt="$a$" src="svgs/44bc9d542a92714cac84e01cbbb7fd61.svg" valign=0.0px width="8.68915409999999pt" height="7.0776222pt"/> equal,
in turn, to say 2, 3, and 5, as would be the case, by Fermat's Little Theorem, if <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> were in fact prime.

Rather than implement Fermat's and related primality tests yourself to detect whether **decimal** is prime,
feel free to use [numlib](https://github.com/sj-simmons/numlib)'s implementation:

```python
import numlib
numlib.isprime(decimal)
```
Exercise 1: Replace the <img alt="$200$" src="svgs/88db9c6bd8c9a0b1527a1cedb8501c55.svg" valign=0.0px width="24.657628049999992pt" height="10.5936072pt"/> above with say <img alt="$k$" src="svgs/63bb9849783d01d91403bc9a5fea12a2.svg" valign=0.0px width="9.075367949999992pt" height="11.4155283pt"/> and write a function that returns a <img alt="$k$" src="svgs/63bb9849783d01d91403bc9a5fea12a2.svg" valign=0.0px width="9.075367949999992pt" height="11.4155283pt"/>-bit prime.

At issue is the fact that larger primes are harder to find. The difficulty is gauged by the prime
number theorem.  If we define <img alt="$\pi(n)$" src="svgs/ab6b1f726144febfe19f0c5d987822fa.svg" valign=-4.109589000000009px width="32.61239849999999pt" height="16.438356pt"/> to be the number of primes less than or equal to <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>, then
the prime number theorem states that <img alt="$\pi(n)$" src="svgs/ab6b1f726144febfe19f0c5d987822fa.svg" valign=-4.109589000000009px width="32.61239849999999pt" height="16.438356pt"/> is well-approximated by <img alt="$n/\ln(n)$" src="svgs/9bbe784fa51a44e5989bc0b2cac489c3.svg" valign=-4.109589000000009px width="57.17672729999999pt" height="16.438356pt"/> in the sense that
<p align="center"><img alt="$$\lim_{n\rightarrow\infty}\pi(n)\cdot\ln(n)/n=1.$$" src="svgs/bbbfdbd8c7959c658b17fe94615640ac.svg" valign=0.0px width="170.7003078pt" height="22.1917806pt"/></p>

The prime number theorem implies that the number of primes between <img alt="$2^k$" src="svgs/91f4e50a1561b60d45e7079ca70f2ed4.svg" valign=0.0px width="15.48523844999999pt" height="13.95621975pt"/> and <img alt="$2^{k+1}$" src="svgs/bf56939689dfdac754c6e27725da93c9.svg" valign=0.0px width="32.12915969999999pt" height="13.95621975pt"/>
is approximately

<p align="center"><img alt="$$\frac{2^{k+1}}{\ln(2^{k+1})}-\frac{2^{k}}{\ln(2^{k})}=\frac{2^k}{\ln(2)}\left(\frac{2}{k+1}-\frac{1}{k}\right)=\frac{2^k}{\ln(2)}\frac{k-1}{k(k+1)};$$" src="svgs/ee03776f623d8b1c2733d4ee1290e882.svg" valign=0.0px width="418.50896339999997pt" height="40.6935375pt"/></p>

hence, the probability of a randomly chosen number between <img alt="$2^k$" src="svgs/91f4e50a1561b60d45e7079ca70f2ed4.svg" valign=0.0px width="15.48523844999999pt" height="13.95621975pt"/> and <img alt="$2^{k+1}$" src="svgs/bf56939689dfdac754c6e27725da93c9.svg" valign=0.0px width="32.12915969999999pt" height="13.95621975pt"/> being
prime is approximately <img alt="$p = (k-1)/(\ln(2)k(k+1))\approx 1/(ln(2)k).$" src="svgs/8a69b5df042e9317a0458345b0d88fc2.svg" valign=-4.109589000000009px width="296.1968427pt" height="16.438356pt"/>
It is an elementary fact from probability theory that, one average, one must test
<img alt="$1/p=\ln(2)k$" src="svgs/d3ae945c3802db1b078977a80a90f615.svg" valign=-4.109589000000009px width="90.40530014999999pt" height="16.438356pt"/> numbers before turning one up that is in fact prime.

Exercise 2: Write a program that verifies the expected number of tries before your
your function from exercise 1 is about 



### Public-key Cryptography

In modern times, you can create and publish (on, say, your personal webpage) a *public key* that
can then be used (by, say, someone called Athena) to encrypt a private message to you.  You can decrypt
Athena's message but no else can, so it doesn't matter if a bad actor sees Athena's encrypted message
that she is sending to you.

Important: since your enciphering key is public, a bad actor might try to intercept Athena's message
and replace it with a malicious message encrypted with your public key.  We need to bar against it
but, for now, let us ignore that weakness.

To determine your public key, you first choose two large primes <img alt="$p$" src="svgs/2ec6e630f199f589a2402fdf3e0289d5.svg" valign=-3.1963502999999895px width="8.270567249999992pt" height="10.2739725pt"/> and <img alt="$q$" src="svgs/d5c18a8ca1894fd3a7d25f242cbe8890.svg" valign=-3.1963502999999895px width="7.928106449999989pt" height="10.2739725pt"/> (which you will keep
secret) and multiply them together.  Your public key is the pair of numbers <img alt="$(e, n)$" src="svgs/1ea602d667e6403959572fffbd4671bc.svg" valign=-4.109589000000009px width="37.61233079999999pt" height="16.438356pt"/> where <img alt="$n=pq$" src="svgs/dd91a3c974e35acb9bfc9b9833a127b8.svg" valign=-3.1963502999999895px width="47.98317974999999pt" height="10.2739725pt"/>
and <img alt="$e$" src="svgs/8cd34385ed61aca950a6b06d09fb50ac.svg" valign=0.0px width="7.654137149999991pt" height="7.0776222pt"/> a positive integer relatively prime to <img alt="$\phi(n)=(p-1)*(q-1)$" src="svgs/5d28570761d98da2134694296c4ca9e6.svg" valign=-4.109589000000009px width="168.27977264999998pt" height="16.438356pt"/>.

Now suppose that <img alt="$M$" src="svgs/fb97d38bcc19230b0acd442e17db879c.svg" valign=0.0px width="17.73973739999999pt" height="11.232861749999998pt"/> is a positive integer representing you version 
If <img alt="$M$" src="svgs/fb97d38bcc19230b0acd442e17db879c.svg" valign=0.0px width="17.73973739999999pt" height="11.232861749999998pt"/> is the numeric version of your message, then we encrypt 

Now, in order for you to decrypt Athena's message you must (see below) derive your
So how is it that your public key can't be reverse engineered by a bad actor

And example

* ECC (Elliptic Curve Cryptography
  * [A gentle intro to ECC](https://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/)
  * Neal Koblitz's 1985 paper [Elliptic Curve Cryptosystems](https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866109-5/S0025-5718-1987-0866109-5.pdf)

* Potentially relevant
  * MC Frontalot's [Secrets of the future](https://www.youtube.com/watch?v=FUPstXCqyus)

#### Cryptocurrencies
* [Chart showing which coins use which flavors of cryptography and curves](http://ethanfast.com/top-crypto.html)

* Literature
  * [The Mathematics of Bitcoin](https://arxiv.org/abs/2003.00001)

#### Zero knowledge
* [math primer](http://extropy.foundation/workshops/zkp/primer.html)

#### More literature
* [evervault.com/papers](https://evervault.com/papers)

## Tools/tutorials

#### Quick start on making basic  number-theoretic/cryptography computations in high-performance Python
* Install the multi-precision library [gmpy2](https://gmpy2.readthedocs.io/en/latest/intro.html)
  in a Debian-like environment such as Ubuntu by issuing, at your commandline, the command
  ```shell
  sudo apt install gmpy2
  ```
* Then, if you say want to generate some large primes in Python:
  ```python
  import gmpy2

  randstate = gmpy2.random_state(1728)

  def getprime(nbits):
      """ Return an n-bit prime of type mpz. """
      return gmpy2.next_prime(gmpy2.mpz_rrandomb(randstate, nbits))
  ```
  The function **getprime()** generates primes. Let us check that they are roughly n-bits.
  ```python
  import math

  for _ in range(10):
     p = getprime(128)
     print(f'{p}  approx. bit-length: {math.log2(p)}')

  # output:
  # 170141183460778792299372374151321878621  approx. bit-length: 127.00000000000263
  # 340281722954356176270385562639172354083  approx. bit-length: 127.9999972697725
  # 340282366920938463463374607431768080437  approx. bit-length: 128.0
  # 340282366920938461102191365996962381853  approx. bit-length: 128.0
  # 340282366919700523719237132310450012183  approx. bit-length: 127.99999999999476
  # 329648544222309736707797093892948492193  approx. bit-length: 127.95419631593471
  # 338953138926372144816642355375684713631  approx. bit-length: 127.99435343686405
  # 340282366920937254539860835811807199379  approx. bit-length: 128.0
  # 170141183460470440657506916146035556433  approx. bit-length: 127.00000000000001
  # 340282366920937292316486855759755210761  approx. bit-length: 128.0
  ```
* Now suppose that we want to (naively) pick two large primes, set up an RSA
  encryption scheme, and encode the message 98765432123456789.
  ```python
  p, q = getprime(128), getprime(128)
  n = p * q  # type mpz

  phi = (p - 1) * (q - 1) # keep this secret

  # find a (naive) encryption key
  e = 1
  while (e < 2 or gmpy2.gcd(e, phi) != 1):
    e = gmpy2.mpz_random(randstate, phi) # type mpz

  assert gmpy2.gcd(e, phi) == 1  # (e, n) is our public key

  m = 98765432123456789 # our message

  c = gmpy2.powmod(m, e, n) # m^e (mod n) our ciphertext (encrypted message)
  ```
* Lastly, let us check things by finding a deciphering key and decrypting the message.
  ```python
  _, d, _ = gmpy2.gcdext(e, phi) # d is our deciphering (private) key
  d = gmpy2.t_mod(d, phi)  # we only need d modulo phi

  # Note: the function gcdext used above applies the extended Euclidean algorithm;
  # i.e., gcdext(e, phi) returns a triple of integers (gcd(d, phi), d, f) satisfying
  #                    gcd(e, phi) = d * e + f * phi.
  # Since gcd(e, phi) = 1, d is just the multiplicative inverse of e (modulo phi)

  assert gmpy2.t_mod(d*e, phi) == 1

  # Alternatively, one can compute d with
  #                    d = gmpy2.inverse(e, phi)

  print(gmpy2.powmod(c, d, n)) # c^d (mod n) = 98765432123456789
  ```
#### Use Simmons' [numlib](https://github.com/sj-simmons/numlib) to instantiate and compute, in Python, in a finite field or in an elliptic curve over a finite field.
*  First install the library by typing, at your commandline,
   ```shell
   pip install polylib --user
   pip install numlib --user
   ```
*  Then, you can create in Python a primefield (i.e., Z/pZ where p is a prime) of order say 2027:
   ```pycon
   >>> import numlib as nl
   >>> PF = nl.Zmodp(2027)
   >>> PF  # Z/2027
   >>> x = PF(100)
   >>> x   # 100 + <2027>
   >>> x**1000  # 2022 + <2027>
   >>> x**-1  # 750 + <2027>
   ```
* Instantiate an work in a Galois field of order 7^3:
  ```pycon
  >>> GF = nl.GaloisField(7, 3)
  >>> GF  # Z/7[t]/<t^3+3t^2-3>
  >>> t = GF.t()
  >>> t   # t + <t^3+3t^2-3>
  >>> t**1000  # -3t^2+3t+2 + <t^3+3t^2-3>
  >>> t**-1  # -2t^2+t + <t^3+3t^2-3>
  ```
* Work in an elliptic curve over the Galois field of order 7^3:
  ```pycon
  >>> E = nl.EllCurve(t**2+5*t+2,t-2)
  >>> E  # y^2 = x^3 + (t^2-2t+2)x + (t-2) over Z/7[t]/<t^3+3t^2-3>
  >>> E.j  # 3t^2+t+3
  >>> # Let us list a few point on this curve:
  >>> finite_points = list(nl.finite(E)) # nl.finite(E) is a generator
  >>> for point in finite_points:
  ...     print(point)
  ...
  (-t^2-t-2, -3t^2-3t-3)
  (-t^2-t-2, 3t^2+3t+3)
  (2t^2-t+3, -2t^2-3t-3)
  (2t^2-t+3, 2t^2+3t+3)
  (3t+2, -3t-3)
       ...
  >>> len(finite_points)  # this curve has order 320 (including the point at infinity)
  319
  >>> # Let's pick point and find its order
  >>> pt = finite_points[0]
  >>> pt  # (-t^2-t-2, -3t^2-3t-3) on y^2 = x^3 + (t^2-2t+2)x + (t-2) over Z/7[t]/<t^3+3t^2-3>
  >>> print(1000 * pt)  # (t^2+2t, 2t^2-3t-2)
  >>> print(-pt)  # (-t^2-t-2, 3t^2+3t+3)
  >>> nl.addorder(pt, exponent = 320)  # 160
  ```

#### Install and use Sage as a Python library
On a debian-like (again, including WSL):
```shell
sudo apt install sagemath
```

#### Install PARI/GP and use it as a Python library
* First, one needs to install [pari-gp](http://pari.math.u-bordeaux.fr/) (so that the
  program **gp**) is in one's PATH.

  On a Debian-like system this is as easy as:
  ```shell
  sudo apt install pari-gp
  ```
  PARI includes its own interactive shell:
  ```shell
  > gp -q
  ?
  ```
  For more, see the PARI/GP [documentation](http://pari.math.u-bordeaux.fr/doc.html) and
  [tutorials](http://pari.math.u-bordeaux.fr/tutorials.html)
* Then there's the python interface to PARI/GP. It's called
  [cypari2](https://github.com/sagemath/cypari2)

## Reference

### Resources
* [Freaking blockchain's: how do they work](https://norswap.com/blockchain-how/)
* [keylength.com](https://www.keylength.com/)

### Libraries
* [Comparison of cryptography libraries](https://en.wikipedia.org/wiki/Comparison_of_cryptography_libraries)

### [Homomophic encryption](https://homomorphicencryption.org/introduction/)
* [HElib's design](https://homenc.github.io/HElib/documentation/Design_Document/HElib-design.pdf)
* [FHEW](https://github.com/lducas/FHEW)
* [Cupcake](https://github.com/facebookresearch/Cupcake) from FacebookResearch; Rust, implements
  lattice-based homomophic encryption
  * [Somewhat practicle homomorphic encryption](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.400.6346&rep=rep1&type=pdf)
* [zama](https://zama.ai/concrete/)
* Google's [FHE](https://github.com/google/fully-homomorphic-encryption)

### Curves
* [Safecurves](https://safecurves.cr.yp.to/index.html)
* [Standard Curve Database](https://neuromancer.sk/std/)

### More
* [Cryptol](https://cryptol.net/) A DSL for cryptography that sits on top of Haskell.
* [Theoretical description and formal security analysis of Signal's protocol](https://eprint.iacr.org/2016/1013)
* [Cryptohack](https://cryptohack.org/)


