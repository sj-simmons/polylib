# polylib

The package provides a class for concretely computing in the ring of polynomials with
coefficients in a given ring or field.  See below for example usage. After installation,
additional documentation is made available by typing in your Python interpreter:
```pycon
>>> import polylib
>>> help(polylib.Polynomial)
>>> help(polylib.FPolynomial)
```
Alternatively, type at your commandline: **pydoc3 polylib.Polynomial**.
For recent versions of Python, pydoc3 is simply pydoc.

Installing this package also installs an executable program, **bernoulli**, that
computes Bernoulli polynomials; type **bernoulli** at your
commandline or, if **~/.local/bin/** is not in your PATH, type **~/.local/bin/bernoulli**.

### Contents
* [Installation](#installation)
* [Basic usage](#basic-usage)
* [Cyclotomic polynomials](#cyclotomic-polynomials)
  * [as polynomials over the integers](#cyclotomic-polynomials-as-polynomials-over-the-integers)
  * [over fields](#cyclotomic-polynomials-over-fields)
  * [over finite fields](#cyclotomic-polynomials-over-finite-fields)
* [Resources](#Resources)
* [Todo](#todo)

## Installation

From your command line:

```shell
pip install polylib --user
```
or, for the latest dev version:
```shell
pip install git+https://github.com/sj-simmons/polylib.git --user
```
For older Python installations, **pip** might refer to its Python2 version
(check this with **pip -V**), so use instead **pip3**.

## Basic usage

In the interpreter, one can
```pycon
>>> import polylib
>>> p = polylib.Polynomial([1, 0, -1])  # define a poly by specifying its coefficients
>>> p  # Polynomial([1, 0, -1])
>>> print(p)  # 1 - x^2
>>> print(p**5)  # 1 - 5x^2 + 10x^4 - 10x^6 + 5x^8 - x^10
```
Equivalently, one can define polynomials in a more
[Sage](https://www.sagemath.org/)-like manner:
```pycon
>>> x = polylib.Polynomial([0, 1])  # an indeterminant from which to build polynomials
>>> p = 1 - x**2
>>> print(p**5)  # 1 - 5x^2 + 10x^4 - 10x^6 + 5x^8 - x^10
```
If you are working over a field, you may wish to use **FPolynomial**
&mdash; the difference is that, in the **divmod**, **floordiv**, and **mod** methods, one
can divide an FPolynomial by any FPolynomial, whereas the divisor must be monic
when dividing **Polynomial**s.
```pycon
>>> from fractions import Fraction as F
>>> from polylib import FPolynomial
>>>
>>> x = FPolynomial([0, F(1)])  # an indeterminant for polys over the rational numbers
>>> p = F(3) - F(7,2)*x - F(3,2)*x**2 - F(11,12)*x**3 - F(5,4)*x**4
>>> print(p)  # 3 - 7/2x - 3/2x^2 - 11/12x^3 - 5/4x^4
>>>
>>> # looking at the constant and leading coefficient of
>>> 12*p  # 36 - 42x - 18x^2 - 11x^3 - 15x^4,
>>> # a possible root of p is 3/5. Let us check:
>>> p.of(F(3,5)) == 0  # True
>>>
>>> # Let us divide p by  x - 3/5  using long division:
>>> p1 = x - F(3,5)
>>> q, r = p.divmod(p1)
>>> print("quotient:", q)   # quotient:  5 + 5/2x + 5/3x^2 + 5/4x^3
>>> print("remainder:", r)  # remainder: 0
>>>
>>> #If you just want the quotient:
>>> print(p // p1)  # 5 + 5/2x + 5/3x^2 + 5/4x^3
>>>
>>> #or, just the remainder:
>>> print(p % p1)  # 0
```
Suppose we want a different name for the indeterminant:
```pycon
>>> from polylib import Polynomial
>>> print(Polynomial([3, 4, 0, 2], 't'))  # 3 + 2t + 4t^3
>>> # equivalently:
>>> t = Polynomial([0,1], 't')
>>> p = 3 + 2*t + 4*t**3  # 3 + 2t + 4t^3
```
#### Polynomials with polynomial coefficients

This library does not robustly handle polynomials in more than one variable; however,
when working, say, over finite fields, we often want polynomials with polynomial
coefficients.

We could do something like this:
```pycon
>>> t = Polynomial([0, 1], 't')
>>> p = Polynomial([5-t, t**2+t-1, -2])
>>> print(p)  # 5 - t + -1 + t + t^2x + -2x^2
```
Problem: the string output is ambiguous; instead, do this:
```pycon
>>> t = Polynomial([0, 1], '(t)')  # note the parentheses
>>> p = Polynomial([5-t, t**2+t-1, -2])
>>> print(p)  # (5 - t) + (-1 + t + t^2)x + -2x^2
```
but that is incorrect since
```pycon
>>> p  # Polynomial((Polynomial((5, -1)), Polynomial((-1, 1, 1)), -2))
```
This is what we want:
```pycon
>>> p = Polynomial([5-t, t**2+t-1, -2*t**0])  # note the use of t**0
>>> print(p)  # (5 - t) + (-1 + t + t^2)x + (-2)x^2
>>> p  # Polynomial((Polynomial((5, -1)), Polynomial((-1, 1, 1)), Polynomial((-2,))))
```
Alternatively, of course, one can type
```pycon
>>> p=Polynomial((Polynomial((5,-1),'(t)'),Polynomial((-1,1,1),'(t)'),Polynomial((-2,),'(t)')))
```
but this might be better:
```pycon
>>> x = Polynomial([0, 1])
>>> t = Polynomial([0, 1], '(t)')  # note the parentheses
>>> p = (5-t) * x**0 + (t**2+t-1) * x - 2*t**0 * x**2  # note, now, both t**0 and x**0
>>> print(p)
(5 - t) + (-1 + t + t^2)x + (-2)x^2
```
These sorts of polynomials occur when working with elliptic curves over a finite field.
Their string representations are more readable if the coefficient polynomials don't appear
with spaces:
```pycon
>>> x = Polynomial([0, 1])
>>> t = Polynomial([0, 1], '(t)', spaces=False)
>>> p = (5-t) * x**0 + (t**2+t-1) * x - 2*t**0 * x**2  # note, now, both t**0 and x**0
>>> print(p)
(5-t) + (-1+t+t^2)x + (-2)x^2
```
In some cases, wrapping in square brackets is more readable:
```pycon
>>> t = Polynomial([0, 1],'[t]')
>>> p = Polynomial([complex(0,1), complex(2,3)-complex(0,1)*t])
>>> print(p)  # 1j + [(2+3j) + (-0-1j)t]x
>>> print(p*2)  # (-1+0j) + [(-6+4j) + (2+0j)t]x + [(-5+12j) + (6-4j)t + (-1+0j)t^2]x^2
```
You may prefer a capital x or some other symbol for the indeterminate:
```pycon
>>> p = Polynomial([complex(0,1), complex(2,3)-complex(0,1)*t], '[X]')
>>> print(p*2)  # (-1+0j) + [(-6+4j) + (2+0j)t]X + [(-5+12j) + (6-4j)t + (-1+0j)t^2]X^2
```

## Cyclotomic polynomials

Next, suppose that we try to compute
[cyclotomic polynomials](https://en.wikipedia.org/wiki/Cyclotomic_polynomial)
by means of their definition in terms of complex
[roots of unity](https://en.wikipedia.org/wiki/Root_of_unity):
```python
from cmath import exp, pi
from math import gcd
from functools import reduce
from operator import mul
from polylib import Polynomial

def cyclotomic_poly(n):
    omega = exp(complex(0,2*pi/n))  # e^(2pi/n i), a primitive n_th root of unity
    x = Polynomial([0, 1])
    return reduce(mul, [x - omega**j for j in range(n) if gcd(j, n) == 1], 1)

print(cyclotomic_poly(8))
```
<a id="complexrootsdef"></a>
The function **cyclotomic_poly(n)** computes, in standard math notation,
<img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/>, which is defined by

<p align="center"><img alt="$$\Phi_n(x) = \prod_{\substack{0&lt;j&lt;n \\ \gcd(j,n)=1}} \left(x - e^{\frac{2\pi j}{n} i}\right),$$" src="svgs/b3e6fbe93ecd3362e655d224b83cd61d.svg" valign=0.0px width="223.1982489pt" height="53.7201786pt"/></p>

the terms in the product on the right involving only *primitive* roots of unity.

The output of the above program is
```shell
(1-1.5543122344752192e-15j) + (7.771561172376096e-16+1.1102230246251565e-16j)x + (3.3306690738754696e-16-2.220446049250313e-16j)x^2 + (6.661338147750939e-16+0j)x^3 + x^4
```

The **cmath** Python library, as well as the builtin type **complex**, represents
complex numbers in rectangular form, as pairs of floats; hence the rounding error
in the coefficients of the 8th cyclotomic polynomial. In fact, as we will see below,
the coefficients of <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> are integral for all <img alt="$n.$" src="svgs/ea8d90fb4a8d92af94283e10af3efb57.svg" valign=0.0px width="14.433101099999991pt" height="7.0776222pt"/>

Since the complex primitive roots of unity occur in complex conjugate pairs.
we can improve things by multiplying first the conjugate pairs:
```python
def cyclotomic_poly(n):
    omega = exp(complex(0,2*pi/n))
    x = Polynomial([0, 1])
    if n < 3:
        poly = reduce(mul, [x - omega**j for j in range(n) if gcd(j, n) == 1], 1)
    else:
        poly = reduce(
            mul,
            [(x-omega**j)*(x-(omega**j).conjugate()) for j in range(n//2+1) if gcd(j, n) == 1],
            1
        )
    return poly
```
Then **cyclotomic_poly(8)** becomes:
```
1 + (-4.440892098500626e-16+0j)x + (4.440892098500626e-16+0j)x^2 + (-4.440892098500626e-16+0j)x^3 + x^4
```
Assuming that the coefficients are integral, we can now simply **round**:

```python
from cmath import exp, pi
from math import gcd
from functools import reduce
from operator import mul
from polylib import Polynomial

def cyclotomic_poly(n):
    omega = exp(complex(0, 2 * pi / n))
    x = Polynomial([0, 1])
    if n < 3:
        poly = reduce(mul, [x - omega ** j for j in range(n) if gcd(j, n) == 1], 1)
    else:
        poly = reduce(
            mul,
            [
                (x - omega ** j) * (x - (omega ** j).conjugate())
                for j in range(n // 2 + 1)
                if gcd(j, n) == 1
            ],
            1,
        )
    return Polynomial([round(coeff.real) for coeff in poly])

print(" n   cyclotomic_poly(n)")
for n in range(1, 17):
    print(f"{n:2}   {cyclotomic_poly(n)}")
```
output:
```shell
 n   cyclotomic_poly(n)
 1   -1 + x
 2   1 + x
 3   1 + x + x^2
 4   1 + x^2
 5   1 + x + x^2 + x^3 + x^4
 6   1 - x + x^2
 7   1 + x + x^2 + x^3 + x^4 + x^5 + x^6
 8   1 + x^4
 9   1 + x^3 + x^6
10   1 - x + x^2 - x^3 + x^4
11   1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10
12   1 - x^2 + x^4
13   1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11 + x^12
14   1 - x + x^2 - x^3 + x^4 - x^5 + x^6
15   1 - x + x^3 - x^4 + x^5 - x^7 + x^8
16   1 + x^8
```
(For those interested in computing *exactly* in the real of complex numbers see,
for example,
[Calcium](https://fredrikj.net/calcium/#python-interface) and
[slides](https://fredrikj.net/math/calcium2021nuscap.pdf) from a talk by
its author.)

In the case of cyclotomic polynomials, it is better to reformulate their
definition completely in terms of <img alt="$\mathbb{Z}[x]$" src="svgs/50222103385d9960679d6dc26ba3c47a.svg" valign=-4.109589000000009px width="29.48637614999999pt" height="16.438356pt"/>, polynomials over the
integers.  Before we do that, let us observe the floating point issues inherent
in the above formulation of **cyclotomic_poly(n)** building up to the
point of throwing off our computations.

It is fairly easy to establish (see below), for <img alt="$p$" src="svgs/2ec6e630f199f589a2402fdf3e0289d5.svg" valign=-3.1963502999999895px width="8.270567249999992pt" height="10.2739725pt"/> a prime and <img alt="$r$" src="svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg" valign=0.0px width="7.87295519999999pt" height="7.0776222pt"/> a positive
integer, the following identity:

<p align="center"><img alt="$$\Phi_{p^r}(x) = \Phi_p\left(x^{p^{r-1}}\right).$$" src="svgs/6a32f49891d60dc3be2c836f25f031ef.svg" valign=0.0px width="157.46044709999998pt" height="29.58934275pt"/></p>

We can use **polylib** to compose polynomials. Let us verify the identity
above for <img alt="$p=3$" src="svgs/39cd1e5f222b0b87f4a26f8b97296bd2.svg" valign=-3.196350299999994px width="38.40740639999999pt" height="13.789957499999998pt"/> and <img alt="$r=4$" src="svgs/0423aef2a8675301835de97021681953.svg" valign=0.0px width="38.009795999999994pt" height="10.5936072pt"/>:
```python
x = Polynomial([0, 1])
cyclotomic_poly(3**4) == cyclotomic_poly(3).of(x**(3**3))  # True
```
But what if <img alt="$r=5$" src="svgs/c4375faa868833f80a54a7668f208f49.svg" valign=0.0px width="38.009795999999994pt" height="10.5936072pt"/>?
```python
cyclotomic_poly(3**5) == cyclotomic_poly(3).of(x**(3**4))  # False, but should be True
```
### Cyclotomic polynomials as polynomials over the integers

Over any field <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/> with unit,
one can construct cyclotomic polynomials <img alt="$\Phi_n(x)\in \mathbb{Z}[x]$" src="svgs/8f9f85baccba44e8114446211d4b335d.svg" valign=-4.109589000000009px width="92.57805314999999pt" height="16.438356pt"/>
(partially) characterized by the following:

<p align="center"><img alt="\begin{align}&#10;\Phi_n(b) = 0 \text{~if and only if~} b \text{~has exponent~} n. \tag{1}&#10;\end{align}" src="svgs/a7778f49308700dedb60031539e6a342.svg" valign=0.0px width="500.8916352pt" height="16.438356pt"/></p>

(Here, <img alt="$b\ne 0$" src="svgs/a38777a8060cb1108625138d18848be2.svg" valign=-3.1963189500000055px width="37.19163689999999pt" height="14.61184725pt"/> and the *exponent* of <img alt="$b\in\mathbb{F}$" src="svgs/92691d245163b943b4ef0c86b033860f.svg" valign=-0.6427030500000053px width="37.19162039999999pt" height="12.05823135pt"/> is the
smallest positive integer <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>, if it exists, satisfying <img alt="$b^n = 1.$" src="svgs/df075450648b7bb6941293bab1493563.svg" valign=0.0px width="50.70579854999998pt" height="11.4155283pt"/>)

More correctly, the coefficients are in an isomorphic image of
<img alt="$\mathbb{Z}$" src="svgs/b9477ea14234215f4d516bad55d011b8.svg" valign=0.0px width="10.95894029999999pt" height="11.324195849999999pt"/> embedded in the prime field <img alt="$\mathbb{Q}$" src="svgs/0f452ec0bcf578fa387e4857f80f03f4.svg" valign=-2.739730950000001px width="12.785434199999989pt" height="14.0639268pt"/> of <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>
in the case that <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/> has characteristic zero; if <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>
has positive characteristic <img alt="$p$" src="svgs/2ec6e630f199f589a2402fdf3e0289d5.svg" valign=-3.1963502999999895px width="8.270567249999992pt" height="10.2739725pt"/>, then the coefficients live in its
prime field &mdash; isomorphic to <img alt="$\mathbb{Z}/p\mathbb{Z}$" src="svgs/a05b826333ec801b65201f7764c6754f.svg" valign=-4.109589000000009px width="38.40765719999999pt" height="16.438356pt"/> &mdash; which
is naturally the homomorphic image of <img alt="$\mathbb{Z}$" src="svgs/b9477ea14234215f4d516bad55d011b8.svg" valign=0.0px width="10.95894029999999pt" height="11.324195849999999pt"/> in <img alt="$\mathbb{F}.$" src="svgs/2b2b9235c734271289d5b75363a19d56.svg" valign=0.0px width="14.61190994999999pt" height="11.324195849999999pt"/>

If <img alt="$\mathbb{F}=\mathbb{C}$" src="svgs/38a89276e20a1ae2efbd72d8763a5491.svg" valign=0.0px width="43.83549554999999pt" height="11.324195849999999pt"/>, notice that, by the
[Fundamental Theorem of Algebra](https://en.wikipedia.org/wiki/Fundamental_theorem_of_algebra),
the cyclotomic polynomials defined [above](#complexrootsdef) satisfy condition (1) since
the primitive <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th roots of unity, <img alt="$e^{\frac{2\pi j}{n}i}$" src="svgs/8c4c460cd4a190cb8a32cf8e829367dc.svg" valign=0.0px width="34.289356199999986pt" height="16.7563605pt"/> for <img alt="$0&lt; j &lt; n$" src="svgs/90e734b49fa3b62a06158cf8953f391d.svg" valign=-3.1963519500000044px width="69.63176549999999pt" height="14.0378568pt"/>
and <img alt="$\gcd(j,n)=1,$" src="svgs/93b5982c2035e80f282812444ef6aa93.svg" valign=-4.109589000000009px width="96.11601449999999pt" height="16.438356pt"/> are the complete list of complex numbers of exponent <img alt="$n.$" src="svgs/ea8d90fb4a8d92af94283e10af3efb57.svg" valign=0.0px width="14.433101099999991pt" height="7.0776222pt"/>
Also note that each primitive <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th root of unity is a zero of mutliplicity 1,
so the <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> in the definition [above](#complexrootsdef) are the unique
smallest-degree, monic polynomials that do the job.

How can we construct the cyclotomic polynomials <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> using only
polynomials of the integers?

First note that, over <img alt="$\mathbb{C}$" src="svgs/81324f07e9ffb7920321df72cc0bee1b.svg" valign=0.0px width="11.87217899999999pt" height="11.324195849999999pt"/>, we have the factorization

<p align="center"><img alt="$$x^n-1=\prod_{0\le j &lt; n} \left(x - e^{\frac{2\pi j}{n} i}\right),$$" src="svgs/b458302bbbdd1c20db250743096676c2.svg" valign=0.0px width="204.17966909999998pt" height="40.8101859pt"/></p>

the product now involving *all* <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th roots of unity &mdash; the factorization
holds since the <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th roots of unity are exactly the <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> roots of <img alt="$x^n-1.$" src="svgs/100caa207d95c5e7ba7d084fdab1423e.svg" valign=-1.3698745500000056px width="51.219550799999986pt" height="12.2895597pt"/>

The <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th roots of unity together comprise a multiplicative group of order <img alt="$n.$" src="svgs/ea8d90fb4a8d92af94283e10af3efb57.svg" valign=0.0px width="14.433101099999991pt" height="7.0776222pt"/>
By Lagrange's Theorem, the order, <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/>, of *any* single <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th root of unity must
divide <img alt="$n.$" src="svgs/ea8d90fb4a8d92af94283e10af3efb57.svg" valign=0.0px width="14.433101099999991pt" height="7.0776222pt"/> Collecting the <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/>th primitive roots of unity for each <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/> that
divides <img alt="$n,$" src="svgs/85bc1f723bdc744666d4f2241b1031f7.svg" valign=-3.1963502999999895px width="14.433101099999991pt" height="10.2739725pt"/> we can write

<p align="center"><img alt="$$x^n-1=(x-1)\prod_{\substack{1 &lt; d \le n \\ d|n}} \Phi_d(x),$$" src="svgs/fe49c898dbbb36b773c023401688b7e8.svg" valign=0.0px width="212.8945533pt" height="52.3497744pt"/></p>

<a id="recursivedef"></a>
where <img alt="$\Phi_d(x)$" src="svgs/afb6c6db132e3039fb05e2d1f71b2c8f.svg" valign=-4.109589000000009px width="41.71759679999999pt" height="16.438356pt"/> are, for now, still defined, as above, in terms of the
primitive <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/>th roots of unity; hence, we can write

<p align="center"><img alt="$$\Phi_n(x)=\frac{x^n-1}{\displaystyle\prod_{\substack{1 \le d &lt; n \\ d|n}} \Phi_d(x)}$$" src="svgs/871cf9ed5410d3d10cf77d093801d7dc.svg" valign=0.0px width="153.41765669999998pt" height="72.58175759999999pt"/></p>

where we have set <img alt="$\Phi_1(x)=x -1.$" src="svgs/0b31a35771b3af5ede73661d6ee5fb1b.svg" valign=-4.109589000000009px width="105.61630364999998pt" height="16.438356pt"/> In fact, we can view this last equality as a
(recursive) *definition* of <img alt="$\Phi_n(x).$" src="svgs/2e6ac9ad1026f2ee7068e2c6e2eb00ec.svg" valign=-4.109589000000009px width="47.56676264999999pt" height="16.438356pt"/>

We now see why <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> has integer coefficients: inductively (with respect to <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>), if we assume that each &mdash; monic, we already know &mdash;  <img alt="$\Phi_d(x)$" src="svgs/afb6c6db132e3039fb05e2d1f71b2c8f.svg" valign=-4.109589000000009px width="41.71759679999999pt" height="16.438356pt"/>
with <img alt="$d|n$" src="svgs/98f9d36d03fb9e804152ab60b8f4d6dc.svg" valign=-4.109589000000009px width="22.989063899999987pt" height="16.438356pt"/> and <img alt="$d&lt;n$" src="svgs/f423eb7544b4a861e759315c0b0a8ef0.svg" valign=-0.6427030500000053px width="40.34047049999999pt" height="12.05823135pt"/> has integral coefficients, the product in the denominator
above is also monic with integer coefficients. The numerator, <img alt="$x^n-1$" src="svgs/0d5d6e639f683b34cdfedc73232a0bbf.svg" valign=-1.3698745500000056px width="46.65332594999999pt" height="12.2895597pt"/>, has integer
coefficients so that long division (we already know that the remainder is 0)
can be executed using only integers.

Rather than writing out the induction proof in detail, let us rewrite
the function **cyclotomic_poly** in a such a way that it uses only polynomials
over the integers:

```python
from functools import reduce
from operator import mul
from polylib import Polynomial

def cyclotomic_poly(n):
    assert n > 0
    x = Polynomial([0, 1])
    return (x ** n - 1) // reduce(
        mul, [cyclotomic_poly(d) for d in range(1, n) if n % d == 0], x ** 0
    )

print(" n   cyclotomic_poly(n)")
for n in range(1, 18):
    print(f"{n:2}   {cyclotomic_poly(n)}")
```
The output is, of course, the same as above,
```shell
 n   cyclotomic_poly(n)
 1   -1 + x
 2   1 + x
 3   1 + x + x^2
 4   1 + x^2
 5   1 + x + x^2 + x^3 + x^4
     ...
```
And the two identities, the second of which we had trouble with above?
```python
x = Polynomial([0, 1])
cyclotomic_poly(3**4) == cyclotomic_poly(3).of(x**(3**3))  # True
cyclotomic_poly(3**5) == cyclotomic_poly(3).of(x**(3**4))  # True
```
Before moving on, let us establish the identity in question, which is

<p align="center"><img alt="$$\Phi_{p^r}(x) = \Phi_p\left(x^{p^{r-1}}\right).$$" src="svgs/6a32f49891d60dc3be2c836f25f031ef.svg" valign=0.0px width="157.46044709999998pt" height="29.58934275pt"/></p>

First notice that, if n is prime, the recursive definition specializes to

<p align="center"><img alt="$$\Phi_p(x)=\frac{x^p-1}{x-1}.$$" src="svgs/cb3c847d12b0606c7db59db717da4ec0.svg" valign=0.0px width="117.3837555pt" height="34.68611685pt"/></p>

More generally,

<p align="center"><img alt="\begin{align}&#10;\Phi_{p^r}(x) &amp;= \frac{x^{p^r}-1}{\Phi_{1}(x)\Phi_{p}(x)\Phi_{p^2}(x)\cdots\Phi_{p^{r-1}}(x)} \notag \\&#10;&amp;= \frac{x^{p^r}-1}{\Phi_{1}(x)\Phi_{p}(x)\Phi_{p^2}(x)\cdots\Phi_{p^{r-2}}(x)\frac{x^{p^{r-1}}-1}{\Phi_{1}(x)\Phi_{p}(x)\Phi_{p^2}(x)\cdots\Phi_{p^{r-2}(x)}}} \notag \\&#10;&amp;= \frac{x^{p^r}-1}{x^{p^{r-1}}-1} \notag \\&#10;&amp;= \Phi_p\left(x^{p^{r-1}}\right). \notag&#10;\end{align}" src="svgs/8a8ebac4e068b993cbe7129ceb575255.svg" valign=0.0px width="467.18886884999995pt" height="183.27877425pt"/></p>

### Cyclotomic polynomials over fields

Over <img alt="$\mathbb{C}$" src="svgs/81324f07e9ffb7920321df72cc0bee1b.svg" valign=0.0px width="11.87217899999999pt" height="11.324195849999999pt"/>, the polynomial <img alt="$x^n-1$" src="svgs/0d5d6e639f683b34cdfedc73232a0bbf.svg" valign=-1.3698745500000056px width="46.65332594999999pt" height="12.2895597pt"/> has <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> distinct roots and so each root
has multiplicity 1; i.e., in its factored form, none of the linear factors
<img alt="$x-\alpha$" src="svgs/ad1baf7900c88bcd43bffd2907b2a45e.svg" valign=-1.3698745499999996px width="40.06268309999999pt" height="10.958925449999999pt"/> of <img alt="$x^n-1$" src="svgs/0d5d6e639f683b34cdfedc73232a0bbf.svg" valign=-1.3698745500000056px width="46.65332594999999pt" height="12.2895597pt"/> are repeated.
It turns out the same is true over any field whose
[characteristic](https://en.wikipedia.org/wiki/Characteristic_(algebra))
does not divide <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> (this follows from the fact that <img alt="$x^n-1$" src="svgs/0d5d6e639f683b34cdfedc73232a0bbf.svg" valign=-1.3698745500000056px width="46.65332594999999pt" height="12.2895597pt"/> and its formal derivative
<img alt="$nx^{n-1}$" src="svgs/bf72f73860980747d361bc4cad8dda72.svg" valign=0.0px width="44.214455999999984pt" height="13.380876299999999pt"/> satisfy <img alt="$1 = n^{-1}x \cdot nx^{n-1}-(x^n-1)$" src="svgs/1e1f23cedaa9debed279fcb962dfe67f.svg" valign=-4.109589000000009px width="203.48543654999997pt" height="17.4904653pt"/> if the characteristic of <img alt="$F$" src="svgs/b8bc815b5e9d5177af01fd4d3d3c2f10.svg" valign=0.0px width="12.85392569999999pt" height="11.232861749999998pt"/> does
not divide <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> and hence are co-prime over <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>)<sup>[1](#footnote1)</sup>.

A natural question to consider is: can we use the
[recursive definition](#recursivedef) from the last section over *any* field?
For completeness, here is the recursive definition:

<p align="center"><img alt="$$\Phi_n(x)=\frac{x^n-1}{\displaystyle\prod_{\substack{1 \le d &lt; n \\ d|n}} \Phi_d(x)}$$" src="svgs/871cf9ed5410d3d10cf77d093801d7dc.svg" valign=0.0px width="153.41765669999998pt" height="72.58175759999999pt"/></p>

At issue is whether we retain the property

<p align="center"><img alt="$$\Phi_n(b) = 0 \text{~if and only if~} b \text{~has exponent~} n$$" src="svgs/4d6743aa6d80c948a3560b77eed403ff.svg" valign=0.0px width="296.94301889999997pt" height="16.438356pt"/></p>

over fields other than the complex numbers.

But the only way that
<img alt="$\alpha\in\mathbb{F}$" src="svgs/6d85e1ff11e20c568a64862f9adeee0f.svg" valign=-0.6427030499999994px width="40.71332264999999pt" height="11.966898899999999pt"/> can satisfy <img alt="$\Phi_n(\alpha) = 0$" src="svgs/dcd6a3f279b64e60f098863a07276969.svg" valign=-4.109589000000009px width="74.31888914999999pt" height="16.438356pt"/> is if <img alt="$\alpha^n=1$" src="svgs/c8b807f3fb4becf7701820495b02b50b.svg" valign=0.0px width="49.66127594999999pt" height="10.91968515pt"/>;
and if the reason <img alt="$\alpha^n=1$" src="svgs/c8b807f3fb4becf7701820495b02b50b.svg" valign=0.0px width="49.66127594999999pt" height="10.91968515pt"/> is that <img alt="$\alpha^d=1$" src="svgs/a79926e57ca1e4a2dae317e601bb0fe9.svg" valign=0.0px width="48.37833329999999pt" height="13.95621975pt"/> for <img alt="$0&lt;d&lt;n$" src="svgs/02df0effbe2654713e6c577b3617238b.svg" valign=-0.6427030500000053px width="70.4773113pt" height="12.05823135pt"/>
(assume without loss of generality that <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/> is the least such positive
integer) then <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/> must divide <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> (by the division algorithm).
But then <img alt="$\alpha$" src="svgs/c745b9b57c145ec5577b82542b2df546.svg" valign=0.0px width="10.57650494999999pt" height="7.0776222pt"/> is a root of the denominator so that <img alt="$x-\alpha$" src="svgs/ad1baf7900c88bcd43bffd2907b2a45e.svg" valign=-1.3698745499999996px width="40.06268309999999pt" height="10.958925449999999pt"/> occurs
as a factor in the denominator. This is a
problem since such a factor in the denominator would cancel with an identical
and non-repeated (if the characteristic of <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/> does not divide <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>) factor
in the numerator.  Hence no such <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/> exists, and the property above holds
over any field.

A subtlety here is that we are now officially working over the prime field of
<img alt="$\mathbb{F}.$" src="svgs/2b2b9235c734271289d5b75363a19d56.svg" valign=0.0px width="14.61190994999999pt" height="11.324195849999999pt"/> Let <img alt="$\overline{\Phi}_n(x)$" src="svgs/6240253b82663cf34884a4318f8ab1eb.svg" valign=-4.109588999999997px width="43.00051634999999pt" height="17.97242865pt"/> denote the cyclotomic polynomial
<img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> with coefficients now in the image of the natural ring homomorphism
<img alt="$\mathbb{Z}\rightarrow\mathbb{F}$" src="svgs/a28e6e67d4f49eccfd3014707fc25a54.svg" valign=0.0px width="46.575228149999994pt" height="11.324195849999999pt"/> (which is reduction modulo
a prime if the characteristic of <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/> is positive and essentially the
identity in characteristic zero).  Clearly, <img alt="$\overline{\Phi}_n(x)=0$" src="svgs/b6d85fb9fd09edab5721d8e32c836029.svg" valign=-4.109588999999997px width="73.13735714999999pt" height="17.97242865pt"/> for
<img alt="$\alpha \in \mathbb{F}$" src="svgs/a0ecc1c428c2f73e66078a059b9e7671.svg" valign=-0.6427030499999994px width="40.71332264999999pt" height="11.966898899999999pt"/> if and only if <img alt="$\Phi_n(\alpha)=0.$" src="svgs/60a24c05e1ee5b4d72ff080295c1489b.svg" valign=-4.109589000000009px width="78.88511234999999pt" height="16.438356pt"/> Using the Euclidean
algorithm, we can write, for <img alt="$\Phi_d(x)$" src="svgs/afb6c6db132e3039fb05e2d1f71b2c8f.svg" valign=-4.109589000000009px width="41.71759679999999pt" height="16.438356pt"/> or any <img alt="$f(x)\in\mathbb{F}[x],$" src="svgs/8cb4f9e7cc9ff9a12a5c6df141ddad35.svg" valign=-4.109589000000009px width="85.22831789999998pt" height="16.438356pt"/>
<img alt="$f(x)=g(x)(x-\alpha) + r$" src="svgs/ba7c6ace5f09e43df0344c805958357d.svg" valign=-4.109589000000009px width="165.3385008pt" height="16.438356pt"/> where <img alt="$r\in \mathbb{F};$" src="svgs/91d3bf74d021f1f0b7f0b8937eb64910.svg" valign=-3.1963502999999998px width="42.57600269999999pt" height="14.52054615pt"/> in fact,
<img alt="$r=f(\alpha)$" src="svgs/54b12c1c3fefb810232180b232d28ac4.svg" valign=-4.109589000000009px width="62.96993009999999pt" height="16.438356pt"/> so that <img alt="$\alpha\in\mathbb{F}$" src="svgs/6d85e1ff11e20c568a64862f9adeee0f.svg" valign=-0.6427030499999994px width="40.71332264999999pt" height="11.966898899999999pt"/> is a root if and only if <img alt="$(x-\alpha)$" src="svgs/484f509c035ee9b410949963858b6f98.svg" valign=-4.109589000000009px width="52.84811069999999pt" height="16.438356pt"/>
divides <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> (which, note, also implies that <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> has at most <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> roots in
<img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>).<sup>[2](#footnote2)</sup>

A useful summary of the preceding argument is that if, in the case that
char(<img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>) does not divide <img alt="$n,$" src="svgs/85bc1f723bdc744666d4f2241b1031f7.svg" valign=-3.1963502999999895px width="14.433101099999991pt" height="10.2739725pt"/> we write

<p align="center"><img alt="$$x^n-1 = \displaystyle\prod_{\substack{1 \le d \le n \\ d|n}} \Phi_d(x),$$" src="svgs/af2555734bfb7a2d559997a6126883fb.svg" valign=0.0px width="159.6640716pt" height="52.3497744pt"/></p>

then the cyclotomic polynomials on the right are pairwise coprime as
polynomials in <img alt="$\mathbb{F}[x].$" src="svgs/b49483c92ae5f8f5ce9eae054c06a5c9.svg" valign=-4.109589000000009px width="33.139345799999994pt" height="16.438356pt"/>

If defined using complex <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th roots of unity, then <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> clearly
has degree <img alt="$\phi(n)$" src="svgs/f4bdf2149704f6b9d6d0068d05021138.svg" valign=-4.109589000000009px width="32.44685399999999pt" height="16.438356pt"/> where <img alt="$\phi$" src="svgs/f50853d41be7d55874e952eb0d80c53e.svg" valign=-3.1963503000000055px width="9.794543549999991pt" height="14.611878599999999pt"/> denotes Euler's phi-function (i.e., <img alt="$\phi(n)$" src="svgs/f4bdf2149704f6b9d6d0068d05021138.svg" valign=-4.109589000000009px width="32.44685399999999pt" height="16.438356pt"/>
is defined to be the number of positive integers less than and coprime
to <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>).

But over any field, the degree is the same, as follows inductively
from the recursive definition and the identity

<p align="center"><img alt="$$\sum_{\substack{1 \le d \le n \\ d|n}} \phi(d) = n.$$" src="svgs/d17c0fff2e0a7e4b599f9d0c4bfe4b30.svg" valign=0.0px width="112.29596564999999pt" height="52.3497744pt"/></p>

Let us reformulate the property above:

If char(<img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>) does not divide <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>, then

<p align="center"><img alt="\begin{align}&#10;\Phi_n(b) = 0 \text{~if and only if~} b \text{~has order~} n \text{~in~}\mathbb{F}^* \tag{1'}&#10;\end{align}" src="svgs/7920c47495ad38127d456459a7f8024f.svg" valign=0.0px width="505.03547985pt" height="16.438356pt"/></p>

where <img alt="$\mathbb{F}^*$" src="svgs/ffa565f555807a0d2e9edbee30e9b852.svg" valign=0.0px width="16.78088114999999pt" height="11.324195849999999pt"/> denotes the non-zero elements (so, units) of the
field <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>.

We are now in a good position to observe the following:

<p align="center"><img alt="$$\text{Every finite subgroup of ~} \mathbb{F}^* text{~is cyclic}.$$" src="svgs/d5cffd88e5f4ae4bf227f58f8aec5e0d.svg" valign=0.0px width="295.7068356pt" height="15.337499099999999pt"/></p>

This is true in any field <img alt="$\mathbb{F}.$" src="svgs/2b2b9235c734271289d5b75363a19d56.svg" valign=0.0px width="14.61190994999999pt" height="11.324195849999999pt"/> Why? Suppose given a subgroup of order
<img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> and consider the identity (in <img alt="$\mathbb{F}[x]$" src="svgs/a420c52aa24a502d60aef830b3b45f9f.svg" valign=-4.109589000000009px width="28.57312259999999pt" height="16.438356pt"/>):

<p align="center"><img alt="$$x^n-1 = \displaystyle\prod_{\substack{1 \le d \le n \\ d|n}} \Phi_d(x);$$" src="svgs/805266daaa077cd6b1165dd7d8eea53a.svg" valign=0.0px width="159.6640716pt" height="52.3497744pt"/></p>

The elements of the given group are exactly the solutions of the polynomials
on the left and hence on the right.  From (1'), we know that the roots of
<img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> are the generators of the given group (of which there are
<img alt="$\phi(n),$" src="svgs/8a840194b0424033754aa8cf28e154d6.svg" valign=-4.109589000000009px width="37.01307719999999pt" height="16.438356pt"/> we also know).  More generally, each (also cyclic, of course)
subgroup of the given group must have order <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/> dividing <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> and so has
the <img alt="$\phi(d)$" src="svgs/b557443bb4cb3f2e9270b661c4fcedfa.svg" valign=-4.109589000000009px width="31.13594054999999pt" height="16.438356pt"/> generators consisting exactly of the roots of
<img alt="$\Phi_n(x).$" src="svgs/2e6ac9ad1026f2ee7068e2c6e2eb00ec.svg" valign=-4.109589000000009px width="47.56676264999999pt" height="16.438356pt"/><sup>[3](#footnote3)</sup>

The picture is now complete. Over any field one has the analogue of the
<img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>th complex roots of unity with exactly the same group structure. The
roots of the cycolotomic polynomials are precisely the generators of
of the subgroups of order <img alt="$d$" src="svgs/2103f85b8b1477f430fc407cad462224.svg" valign=0.0px width="8.55596444999999pt" height="11.4155283pt"/> dividing <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>.

In a field of characteristic zero, one has finite groups of all orders
and they are all cyclic. The same is true in finite field, at least for
orders not divisible by the characteristic.

In any case, to find a generator, or *primitive* element,  of a subgroup
of order <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>, one only need to find a root of <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/>.

The program below uses [numlib](https://pypi.org/project/numlib/)
to instantiate the prime field <img alt="$\mathbb{Z}/2\mathbb{Z}$" src="svgs/8c9f3fec52ab01f5d310bfd440cdb199.svg" valign=-4.109589000000009px width="38.35630094999999pt" height="16.438356pt"/> and then
computes the first few cyclotomic polynomials with coefficients
reduced modulo 2.

```python
from functools import reduce
from operator import mul
from polylib import Polynomial
from numlib import Zmod

def cyclotomic_poly(n: int, x: Polynomial = Polynomial([0, 1])) -> Polynomial:
    return (x ** n - 1) // reduce(
        mul, [cyclotomic_poly(d) for d in range(1, n) if n % d == 0], x ** 0
    )

Zp = Zmod(2)
x = Polynomial([0, Zp(1)])

print(" n   n_th cyclotomic poly over", Zp)
for n in range(1, 8):
    print(f"{n:2}   {cyclotomic_poly(n, x)}")
```
output:
```shell
 n   n_th cyclotomic poly over Z/2Z
 1   1 + x
 2   1 + x
 3   1 + x + x^2
 4   1 + x^2
 5   1 + x + x^2 + x^3 + x^4
 6   1 + x + x^2
 7   1 + x + x^2 + x^3 + x^4 + x^5 + x^6
```

### Cyclotomic polynomials over fields


Cyclotomic polynomials are often in the mix when constructing
finite (Galois) fields. Suppose <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> is an irreducible polynomial over
<img alt="$\mathbb{Z}/p\mathbb{Z}$" src="svgs/a05b826333ec801b65201f7764c6754f.svg" valign=-4.109589000000009px width="38.40765719999999pt" height="16.438356pt"/> of degree <img alt="$m.$" src="svgs/4bc868b30aee0dfd5de42ea15b2cb2d8.svg" valign=0.0px width="18.99932429999999pt" height="7.0776222pt"/>  Then,
<img alt="$\mathbb{Z}/p\mathbb{Z} [x]$" src="svgs/faa721d8e506adcff9df5396c42b78a8.svg" valign=-4.109589000000009px width="56.93509304999999pt" height="16.438356pt"/> modulo the ideal <img alt="$\langle f(x)\rangle$" src="svgs/c6084175338a8a0f04a0ff398ca276b0.svg" valign=-4.109589000000009px width="44.78326709999998pt" height="16.438356pt"/> generated
by <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> is a finite field of order <img alt="$p^m$" src="svgs/71b130cb27f4df588e04a4959bb68240.svg" valign=-3.1963519500000044px width="19.935417149999992pt" height="14.116037099999998pt"/>. See, for example,
[numlib](https://github.com/sj-simmons/numlib) for the mathematical
details and implementations.

Let us consider the case in which <img alt="$p=2$" src="svgs/90264925fb137831c8f410cd14c75cff.svg" valign=-3.196350299999994px width="38.40740639999999pt" height="13.789957499999998pt"/> and <img alt="$m=4$" src="svgs/5abf8a13476a371b544bfce3ff081f70.svg" valign=0.0px width="44.56994024999999pt" height="10.5936072pt"/>. Using **numlib** we can
construct the Galois field of order 16 using **numlib.FPmod** but
we need a quartic that is irreducible over the prime field 
<img alt="$\mathbb{Z}/2\mathbb{Z}$" src="svgs/8c9f3fec52ab01f5d310bfd440cdb199.svg" valign=-4.109589000000009px width="38.35630094999999pt" height="16.438356pt"/> (more on finding such a polynomial below).


```pycon
>>> from numlib import Zmod, FPmod
>>> from polylib import FPolynomial
>>>
>>> PF = Zmod(2)  # the prime field
>>> t = FPolynomial([0, PF(1)], 't')  # an indeterminant for Z/2Z[x]
>>> irred = 1 + t**3 + t**4  # this is irreducible over Z/2Z[x]
>>> GF = FPmod(irred)  # representation of the Galois field of order 16
>>> print(GF)  # Z/2Z [t] / <1 + t^3 + t^4>
```
We already know that the multiplicative group of units in **GF** is
cyclic, and that its generators are exactly the roots of the cyclotomic
polynomial <img alt="$\Phi_{15}$" src="svgs/af89a6736f2fbbe7e580da966db5191e.svg" valign=-2.4657286499999893px width="24.97727264999999pt" height="13.698590399999999pt"/>.
As a sanity check, let us brute force verify this. Continuing the previous
interactive session:
```pycon
>>> def order(elt):   # return the multiplicative orde
...     for j in range(1,16):
...         if elt**j == 1:
...             return j
>>>
>>> from itertools import product
>>>
>>> generators = []
>>> for coeffs in product(PF, repeat=4):
...     elt = GF(coeffs)
...     if order(elt) == 15:
...         generators.append(elt)
>>> print(', '.join([str(gen).replace(" ","") for gen in generators]))
t^2, t^2+t^3, t, t+t^2, t+t^2+t^3, 1+t^3, 1+t^2+t^3, 1+t+t^2
```
Now for the verification:
```pycon
>>> from polylib import cyclotomic
>>> cyclo15 = cyclotomic(15, t)  # 1 + t + t^3 + t^4 + t^5 + t^7 + t^8
>>> all(map(lambda elt: cyclo15.of(elt) == 0, generators))
True
```
Serendipitously, the choice of <img alt="$1+t^3+t^4$" src="svgs/2bfc9d2d97401d5491e2603d9d557049.svg" valign=-1.3698729000000083px width="74.20079369999999pt" height="14.750749199999998pt"/> for the irreducible when
implementing <img alt="$\operatorname{GF}(2^4)$" src="svgs/8f8e33b4f9c783ef5910665a99d66f0b.svg" valign=-4.109589000000009px width="52.00929029999999pt" height="17.4904653pt"/> led to the equivalence class
<img alt="$t + \langle 1+t^3+t^4\rangle$" src="svgs/69eb5cd6c9b9a986343ddd29198ffc0c.svg" valign=-4.109589000000009px width="113.83542885pt" height="17.4904653pt"/> being a generator.  That happening
is particularly convenient and we say that <img alt="$1+t^3+t^4$" src="svgs/2bfc9d2d97401d5491e2603d9d557049.svg" valign=-1.3698729000000083px width="74.20079369999999pt" height="14.750749199999998pt"/> is a
*primitive polynomial* for <img alt="$\operatorname{GF}(2^4)$" src="svgs/8f8e33b4f9c783ef5910665a99d66f0b.svg" valign=-4.109589000000009px width="52.00929029999999pt" height="17.4904653pt"/>; i.e.,
an irreducible <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> is a primitive polynomial for
<img alt="$\operatorname{GF}(2^4)$" src="svgs/8f8e33b4f9c783ef5910665a99d66f0b.svg" valign=-4.109589000000009px width="52.00929029999999pt" height="17.4904653pt"/> if <img alt="$t$" src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg" valign=0.0px width="5.936097749999991pt" height="10.110901349999999pt"/>'s equivalence class is a generator
for <img alt="$\mathbb{Z}/2\mathbb{Z}[x]/\langle 1+t^3+t^4 \rangle$" src="svgs/2463781e8bfdf91dd19ea58a46e18ede.svg" valign=-4.109589000000009px width="152.91108524999999pt" height="17.4904653pt"/>.

We know that we need an quartic irreducible over <img alt="$\mathbb{Z}/2\mathbb{Z}[x]$" src="svgs/14f910bc51b9d1d442f82027720000e1.svg" valign=-4.109589000000009px width="56.88373679999999pt" height="16.438356pt"/>
to even instantiate <img alt="$\operatorname{GF}(2^4)$" src="svgs/8f8e33b4f9c783ef5910665a99d66f0b.svg" valign=-4.109589000000009px width="52.00929029999999pt" height="17.4904653pt"/>; but how can we choose one
so that <img alt="$t$" src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg" valign=0.0px width="5.936097749999991pt" height="10.110901349999999pt"/> turns out to be a generator?  A hint is that the irreducible
we chose above divides <img alt="$\Phi_{15}$" src="svgs/af89a6736f2fbbe7e580da966db5191e.svg" valign=-2.4657286499999893px width="24.97727264999999pt" height="13.698590399999999pt"/>:

```pycon
>>> print(cyclo15 % irred)
0
>>> print(cyclo15 // irred)
1 + t + t^4
```
Moreover, <img alt="$1 + t + t^4$" src="svgs/360eb97978ab7370ce4d7d31fd005a5e.svg" valign=-1.3698729000000083px width="66.82633319999998pt" height="14.750749199999998pt"/> is also irreducible in <img alt="$\mathbb{Z}/2\mathbb{Z}[x]$" src="svgs/14f910bc51b9d1d442f82027720000e1.svg" valign=-4.109589000000009px width="56.88373679999999pt" height="16.438356pt"/>.

But how do we find a candidate for a primitive polynomial in the
first place? One clue is that, in the previous example, <img alt="$\Phi_{15}$" src="svgs/af89a6736f2fbbe7e580da966db5191e.svg" valign=-2.4657286499999893px width="24.97727264999999pt" height="13.698590399999999pt"/>
factors like

begin(align)
1+t+t^3+t^4+t^5+t^7+t^8 &= 1 + t^3 + t^4 + t + t^4(t + t^3 + t^4) \notag\\
&= 1 + (1 + t^4)(t + t^3 + t^4) \notag\\

### Exercises

1. Write a program that find an primitive polynomial for
   $\operatorname{GF}(2^12).$




## Footnotes
<a name="footnote1"><sup>1</sup></a> The *formal* (or *algebraic*) derivative
<img alt="$D: \mathbb{F}[x] \rightarrow \mathbb{F}[x],$" src="svgs/5b222387c4fe16da0538c032cd687265.svg" valign=-4.109589000000009px width="115.04772179999998pt" height="16.438356pt"/> which one first defines on monomials as
<img alt="$D(x^n)=nx^{n-1}$" src="svgs/94ec1e856723816fb8bac34f20558a9d.svg" valign=-4.109589000000009px width="111.32667644999998pt" height="17.4904653pt"/> and then extends linearly, satisfies the familiar product rule
from Calculus.  Hence, if <img alt="$f(x) = (x-\alpha)^2\cdot g(x),$" src="svgs/44935f6d43e92c803f80226595cf00c2.svg" valign=-4.109589000000009px width="161.18701994999998pt" height="17.4904653pt"/> then
<img alt="$D[f(x)]=(x-\alpha)^2\cdot D[g(x)] + 2(x-\alpha)\cdot g(x)$" src="svgs/a540ede6a14a9f643d47ee762caa2cf2.svg" valign=-4.109589000000009px width="326.6594265pt" height="17.4904653pt"/>; so that <img alt="$(x-\alpha)$" src="svgs/484f509c035ee9b410949963858b6f98.svg" valign=-4.109589000000009px width="52.84811069999999pt" height="16.438356pt"/>
is a common divisor of <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> of <img alt="$D[f(x)].$" src="svgs/f00e697d424277ff60b17ddd8ef6355f.svg" valign=-4.109589000000009px width="59.76273764999999pt" height="16.438356pt"/>

Also, that the expression <img alt="$1 = n^{-1}x \cdot nx^{n-1}-(x^n-1)$" src="svgs/1e1f23cedaa9debed279fcb962dfe67f.svg" valign=-4.109589000000009px width="203.48543654999997pt" height="17.4904653pt"/> implies that
<img alt="$nx^{n-1}$" src="svgs/bf72f73860980747d361bc4cad8dda72.svg" valign=0.0px width="44.214455999999984pt" height="13.380876299999999pt"/> and <img alt="$x^n-1$" src="svgs/0d5d6e639f683b34cdfedc73232a0bbf.svg" valign=-1.3698745500000056px width="46.65332594999999pt" height="12.2895597pt"/> have no common factors tacitly uses the fact that
<img alt="$\mathbb{F}[x]$" src="svgs/a420c52aa24a502d60aef830b3b45f9f.svg" valign=-4.109589000000009px width="28.57312259999999pt" height="16.438356pt"/> is a
[unique factorization domain](https://en.wikipedia.org/wiki/Unique_factorization_domain)
&mdash; more computationally, <img alt="$\mathbb{F}[x]$" src="svgs/a420c52aa24a502d60aef830b3b45f9f.svg" valign=-4.109589000000009px width="28.57312259999999pt" height="16.438356pt"/> is a
[Euclidean domain](https://en.wikipedia.org/wiki/Euclidean_domain) and,  in fact,
one in which polynomial division can be used to implement the Euclidean algorithm.

<a name="footnote2"><sup>2</sup></a> The ring <img alt="$\mathbb{Z}[x]$" src="svgs/50222103385d9960679d6dc26ba3c47a.svg" valign=-4.109589000000009px width="29.48637614999999pt" height="16.438356pt"/>, where the cyclotomic
polynomials <img alt="$\Phi_n(x)$" src="svgs/e6acced1a4665534fce181060a1a7b7c.svg" valign=-4.109589000000009px width="43.00053779999999pt" height="16.438356pt"/> live, is not a Euclidean domain, nor even an principal ideal domain;
it is, however, a unique factorization domain in which the division algorithm can be
used to write any <img alt="$f(x)$" src="svgs/7997339883ac20f551e7f35efff0a2b9.svg" valign=-4.109589000000009px width="31.99783454999999pt" height="16.438356pt"/> in <img alt="$\mathbb{Z}[x]$" src="svgs/50222103385d9960679d6dc26ba3c47a.svg" valign=-4.109589000000009px width="29.48637614999999pt" height="16.438356pt"/> as <img alt="$f(x)=g(x)(x-a) + r$" src="svgs/1b1ed27418bdfd527e9ed9ebe217cd84.svg" valign=-4.109589000000009px width="163.45115654999998pt" height="16.438356pt"/> where <img alt="$r\in \mathbb{Z},$" src="svgs/179b9a85f34653eab0cc78fce58dc665.svg" valign=-3.1963502999999998px width="43.489257899999984pt" height="14.52054615pt"/>
so that <img alt="$a\in \mathbb{Z}$" src="svgs/44fab3c024d4b09c2321e3290c56b6cd.svg" valign=-0.6427030499999994px width="39.73923194999999pt" height="11.966898899999999pt"/> is a root of <img alt="$f$" src="svgs/190083ef7a1625fbc75f243cffb9c96d.svg" valign=-3.1963503000000055px width="9.81741584999999pt" height="14.611878599999999pt"/> if and only if it is divisible by <img alt="$x-a$" src="svgs/95e9d73508dbafa46224fc9733ebf0d5.svg" valign=-1.3698745499999996px width="38.17533389999999pt" height="10.958925449999999pt"/>.

<a name="footnote3"><sup>3</sup></a>
If <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/> is [perfect](https://encyclopediaofmath.org/wiki/Perfect_field)
and char(<img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>) divides <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/>, then one can
show that <img alt="$x^n-1$" src="svgs/0d5d6e639f683b34cdfedc73232a0bbf.svg" valign=-1.3698745500000056px width="46.65332594999999pt" height="12.2895597pt"/> necessarily has repeated roots in <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/>. In other words,
there are not enough element of <img alt="$\mathbb{F}$" src="svgs/2d4c6ac334688c42fb4089749e372345.svg" valign=0.0px width="10.045686749999991pt" height="11.324195849999999pt"/> with the right order to form
a group of order <img alt="$n$" src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" valign=0.0px width="9.86687624999999pt" height="7.0776222pt"/> &mdash; so over, for example, over a finite field, there
is no subgroup with order divisible by the characteristic.

## Resources
* Victor Shoup's [A Computational Introduction to Number Theory and Algebra](https://shoup.net/ntb/)
* [FLINT](http://flintlib.org/development.html) Sage uses this for polynomial arithmetic factoring over <img alt="$\mathbb{Z}$" src="svgs/b9477ea14234215f4d516bad55d011b8.svg" valign=0.0px width="10.95894029999999pt" height="11.324195849999999pt"/>, <img alt="$\mathbb{Q}$" src="svgs/0f452ec0bcf578fa387e4857f80f03f4.svg" valign=-2.739730950000001px width="12.785434199999989pt" height="14.0639268pt"/>, and <img alt="$\mathbb{Z}/n\mathbb{Z}$" src="svgs/94d333ba0aaa5e9c8ce88690986075c2.svg" valign=-4.109589000000009px width="40.00396784999999pt" height="16.438356pt"/>.
* Arnold and Monagan's article [Calculating Cyclotomic Polynomials](https://www.ams.org/journals/mcom/2011-80-276/S0025-5718-2011-02467-1/S0025-5718-2011-02467-1.pdf)
* Brent's [On computing factors of cyclotomic polynomials](https://maths-people.anu.edu.au/~brent/pd/rpb135.pdf)
* Cao and Cao's [Note on fast division algorithm for polynomials using Newton iteration](https://arxiv.org/pdf/1112.4014.pdf)

## Todo
* Implement faster division algorithm (see Cao and Cao's paper linked to above).
* Implement multivariate polynomials (use hash?)
* Add option to **Polynomial**'s constructor to print in decreasing degree order.


