#!/usr/bin/env python
"""
Computes the nth Bernoulli number by computing the requisite number of
terms of the generating series x/(1-e^(-x)). (This is not an efficient
way to compute Bernoulli numbers.)                            -Simmons

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

  py -m timeit -s "import bernoulli" "bernoulli.berni(80)"
"""

__author__ = 'Scott Simmons'
__version__ = '0.1'
__status__ = 'Development'
__date__ = '5/1/21'
__copyright__ = """
  Copyright 2014-2021 Scott Simmons

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
__license__= 'Apache 2.0'

import sys
import math
from fractions import Fraction
from polylib import Polynomial
from numbers import Integral


def berniPoly(n):
    """Return the generating series x/(1-e^(-x)) modulo (x^(n+1))
       as a Polynomial.

       One has: Sum_{i=0} B_n (-x)^n/n! = x/(1-e^(-x)).

    Interactively, e.g.,

    >>> print(berniPoly(10))
    1 + 1/2x + 1/12x^2 - 1/720x^4 + 1/30240x^6 - 1/1209600x^8 + 1/47900160x^10

    The corresponding Bernoulli numbers:

    >>> for i,x in enumerate(berniPoly(10)):
    ...    print("B_"+str(i)+" =",(-1)**i*math.factorial(i)*x,end=', ')
    ...                                                    # doctest: +ELLIPSIS
    B_0 = 1, B_1 = -1/2, B_2 = 1/6, ...

    """
    # generate n terms of the polynomial p(x) satisfying x/(1-e^(-x)) = 1/(1-p(x))
    p = [0]
    for i in range(2,n+2):
        p.append(Fraction((-1)**i,math.factorial(i)))
    p = Polynomial(p)

    q = Polynomial([Fraction(1)])

    for i in range(1,n+1):
        q = q *p + Polynomial([1])
        q = Polynomial(q[:n+2])

    return Polynomial(q[:n+1])


def berni(n):
    """Return B_n, the nth Bernoulli number.

       B_n is defined by:

             Sum_{i=0} B_n (-x)^n/n! = x/(1-e^(-x)).

    Interactively, e.g.,

    >>> print(berni(16))
    -3617/510

    """
    q = berniPoly(n)

    # be careful: Polynomial strips trailing zeros
    if q.degree() < n:
        return 0
    else:
        return (-1)**n*math.factorial(n)*q[n]

def main():

    if len(sys.argv) == 2 and sys.argv[1] == '-v':
        import doctest
        doctest.testmod(verbose=True)
        sys.exit()

    if not ( 2 <= len(sys.argv) <= 3 ):
        sys.exit(print(__doc__))

    n = None; show = False

    for arg in sys.argv:    # process command line
        if 's' in arg:
            show = True
        if arg.isdigit():
            n = int(arg)

    if n is None or not isinstance(n,Integral) or n < 0:
        sys.exit(print(__doc__))

    p = berniPoly(n)

    if show == True:
        print("\nx/(1-e^(-x)) =\n",p,"+ O(x^"+str(n+1)+")")
        if p.degree() < n:
            print("\nB_"+str(n)," = ",0)
        else:
            print("\nB_"+str(n)," = ",(-1)**n*math.factorial(n)*p[-1])

    else:
        print("\nB_"+str(n)," = ",berni(n))

if __name__ == '__main__':

    main()
