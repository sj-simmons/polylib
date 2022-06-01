#!/usr/bin/env python
"""
Computes the nth Bernoulli number by computing the requisite number of
terms of the generating series x/(1-e^(-x)). (This is not an efficient
way to compute Bernoulli numbers.)

  Usage: py bernoulli.py [options] n

    arguments:

        n  non-negative integer

    options:

       -s  show x/(1-e^(x)) modulo x^n

       -v  if called with only this option, run tests and exit

  or, interactively, e.g.:

    >>> from bernoulli import *
    >>> print(berni(12))
    -691/2730

    >>> print(berniPoly(8))
    1 + 1/2x + 1/12x^2 - 1/720x^4 + 1/30240x^6 - 1/1209600x^8

Note: Simmons uses this to gauge performance of polynomial computations
for various implementations. For timing, try something like:

  python -m timeit -s "from polylib.bernoulli import berni" "berni(80)"
"""

__author__ = "Scott Simmons"
__version__ = "0.2"
__status__ = "Development"
__date__ = "03/25/22"
__copyright__ = """
  Copyright 2014-2022 Scott Simmons

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
"""
__license__ = "Apache 2.0"

import sys
import math
from fractions import Fraction
from polylib.polynomial import Polynomial


def berniPoly(n: int) -> Polynomial[Fraction]:
    """Return the generating series x/(1-e^(-x)) modulo (x^(n+1))
       as a Polynomial.

       One has: Sum_{i=0} B_n (-x)^n/n! = x/(1-e^(-x)).

    Interactively, e.g.,

    >>> print(berniPoly(10))
    1 + 1/2x + 1/12x^2 - 1/720x^4 + 1/30240x^6 - 1/1209600x^8 + 1/47900160x^10

    The corresponding Bernoulli numbers:

    >>> for i,x in enumerate(berniPoly(10)):
    ...    print("B_"+str(i)+" =",(-1)**i*math.factorial(i)*x,end=', ')
    ...    # doctest: +ELLIPSIS
    B_0 = 1, B_1 = -1/2, B_2 = 1/6, ...

    """
    # generate n terms of the series p(x) satisfying x/(1-e^(-x)) = 1/(1-p(x))
    p_ = [Fraction(0)]
    for i in range(2, n + 2):
        p_.append(Fraction((-1) ** i, math.factorial(i)))
    p = Polynomial(p_)

    # q = Polynomial([Fraction(1)])
    # for i in range(1, n + 1):
    #    q = q * p + Polynomial([Fraction(1)])
    #    q = Polynomial(q._coeffs[: n + 2])

    # return Polynomial(q._coeffs[: n + 1])

    return Polynomial((Fraction(1),)) if n == 0 else (1 - p).formalinv(n)


def berni(n: int) -> Fraction:
    """Return B_n, the nth Bernoulli number.

       B_n is defined by:

             Sum_{i=0} B_n (-x)^n/n! = x/(1-e^(-x)).

    Interactively, e.g.,

    >>> print(berni(16))
    -3617/510

    """
    assert n >= 0
    q = berniPoly(n)

    # be careful: Polynomial strips trailing zeros
    if q._degree < n:
        return Fraction(0)
    else:
        return Fraction((-1) ** n * math.factorial(n) * q[n])


def main() -> None:

    if len(sys.argv) == 2 and sys.argv[1] == "-v":
        import doctest

        # doctest.testmod()
        doctest.testmod(verbose=False)
        sys.exit()

    if not (2 <= len(sys.argv) <= 3):
        sys.exit(__doc__)

    n = None
    show = False

    for arg in sys.argv[1:]:  # process command line
        if "s" in arg:
            show = True
        if arg.isdigit():
            n = int(arg)

    if n is None or not isinstance(n, int) or n < 0:
        sys.exit(__doc__)

    p = berniPoly(n)

    if show == True:
        print("\nx/(1-e^(-x)) =\n", p, f"+ O({str(Polynomial([0]*(n+1)+[1]))})")
        if p._degree < n:
            print("\nB_" + str(n), " = ", 0)
        else:
            print("\nB_" + str(n), " = ", (-1) ** n * math.factorial(n) * p[-1])

    else:
        print("B_" + str(n), " = ", berni(n))


if __name__ == "__main__":

    main()
