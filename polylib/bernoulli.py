#!/usr/bin/env python
"""
Computes the nth Bernoulli number by computing the requisite number of
terms of the generating series x/(1-e^(-x)). This is not an efficient
way to compute Bernoulli numbers. It computes all Bernoulli numbers up
to the nth one.

  Usage: py bernoulli.py [options] n

    arguments:

        n   A non-negative integer.

    options:

       -s   Show x/(1-e^(x)) modulo x^(n+1).

  or, interactively, e.g.:

    >>> from polylib.bernoulli import berni_num
    >>> print(berni_num(12))
    -691/2730
    >>> print(berni_num(40))
    -261082718496449122051/13530

  or:

    >>> from polylib.polynomials import bernoulli
    >>> print(bernoulli(8))
    1 + 1/2x + 1/12x^2 - 1/720x^4 + 1/30240x^6 - 1/1209600x^8

Note: Simmons uses this to gauge performance of various rational number
implementations. For timing, try something like:

python -m timeit -s "from polylib.bernoulli import berni_num as b" "b(500)"
"""

__author__ = "Scott Simmons"
__version__ = "0.3.1"
__status__ = "Development"
__date__ = "01/27/24"
__copyright__ = """
  Copyright 2014-2024 Scott Simmons

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
from polylib.polynomials import Rational, bernoulli


def berni_num(n: int) -> Rational:
    """
    Return B_n, the nth Bernoulli number.

    B_n is defined by: Sum_{i=0} B_n (-x)^n/n! = x/(1-e^(-x)).

    This computes the generating series to the requisite term and
    is therefore not an efficient way to Bernoulli numbers.

    Interactively,

    >>> print(berni_num(16))
    -3617/510
    """
    assert n >= 0
    q = bernoulli(n)

    if len(q.cos) <= n:
        return Rational(0)
    else:
        return Rational((-1) ** n * math.factorial(n)) * q.cos[n]


def main() -> None:

    if not (2 <= len(sys.argv) <= 3):
        sys.exit(__doc__)

    n = None
    show = False

    for arg in sys.argv[1:]:
        if "s" in arg:
            show = True
        if arg.isdigit():
            n = int(arg)

    if n is None or not isinstance(n, int) or n < 0:
        sys.exit(__doc__)

    p = bernoulli(n)  # note: this is computed even if n is odd

    if show == True:
        #if p._degree < n:
        if len(p.cos) <= n:
            print(f"B_{str(n)} = 0")
        else:
            print(f"B_{str(n)} = {Rational((-1) ** n * math.factorial(n)) * p.cos[-1]}")
        print(f"x/(1-e^(-x)) = {p} + O(x^{n+1})")
    else:
        print(f"B_{str(n)} = {berni_num(n)}")


if __name__ == "__main__":

    if len(sys.argv) == 2 and sys.argv[1] == "--doctest":
        import doctest

        doctest.testmod(verbose=False)
        sys.exit()

    main()
