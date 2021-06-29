from polylib import Polynomial
from functools import reduce
from operator import mul

__author__ = "Scott Simmons"
__version__ = "0.1"
__status__ = "Development"
__date__ = "06/23/21"
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
__license__ = "Apache 2.0"


def cyclotomic(n: int, x: Polynomial = Polynomial([0, 1])) -> Polynomial:
    """Return the nth cyclotomic polynomial.

    For a positive integer n, this returns the nth cyclotomic polynomial with,
    by default int coefficients.

    Examples:

        >>> print(cyclotomic(1))
        -1 + x
        >>> print(cyclotomic(2))
        1 + x
        >>> print(repr(cyclotomic(2)))
        Polynomial((1, 1))

        If you want int coefficients but of type Fraction:

        >>> from polylib import FPolynomial
        >>> from fractions import Fraction
        >>> x = FPolynomial([0, Fraction(1)]) # an indeterminant for Q[x]
        >>> print(cyclotomic(2, x))
        1 + x
        >>> print(repr(cyclotomic(2, x)))
        FPolynomial((Fraction(1, 1), Fraction(1, 1)))

        If you want to reduce the coefficients modulo n:

        >>> from polylib import FPolynomial
        >>> from numlib import Zmod
        >>> GF = Zmod(17) # Z/17Z
        >>> x = FPolynomial([0, GF(1)]) #  indeterminant for Z/17Z[x]
        >>> print(cyclotomic(2, x))
        1 + x
        >>> print(repr(cyclotomic(2, x)))
        FPolynomial((1 + <17>, 1 + <17>))
    """
    if not isinstance(n, int):
        raise TypeError(f"n must be of type int, not {type(n)}")
    if n < 1:
        raise ValueError(f"n must be positive, not {n}")

    return (x ** n - 1 * x ** 0) // reduce(
        mul, [cyclotomic(d) for d in range(1, n) if n % d == 0], x ** 0
    )


if __name__ == "__main__":

    import doctest

    doctest.testmod()
