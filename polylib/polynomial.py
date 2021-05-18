"""Provide polynomial classes.

 The base class here is a polynomial class implemented using composition with tuple objects.
 The underlying concrete representation of a polynomial is a tuple of its coefficients of
 length equal to its degree.

 TODO:

 - When we can assume Python10, no need to import annotations from __future__ (which is cur-
   rently used for the type hints for other).
 - If/when Python's Number numeric hierarchy works with type hints, one might be able to use
   nominal typing instead of structural typing.
   Or, or might be better to keep structural typing in any case (so we can use gmpy2 or numpy
   ints or any numbers that cannot be made to subclass Number).
 - Later verions os mypy import Protocal from typing, not typing_extensions
"""

from __future__ import annotations
import copy
from numbers import Number
from typing import Sequence, TypeVar, Union, Tuple, List, Optional
from typing_extensions import Protocol, runtime_checkable

__author__ = "Scott Simmons"
__version__ = "0.1"
__status__ = "Development"
__date__ = "5/1/21"
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

T = TypeVar('T')
#S = TypeVar('S', covariant = True)

@runtime_checkable
class Ring(Protocol):

    """The minimal operations that we require for a ring.

    This in only used for structural typing purposes.
    """

    def __add__(self: T, other: T) -> T: ...

    def __neg__(self: T) -> T: ...

    def __mul__(self: T, other: T) -> T: ...

    def __pow__(self: T, n: int) -> T: ...

    #def __eq__(self: T, other: object) -> bool: ...

class Polynomial:

    """Implements addition, subtraction, multiplication, and evaluation of polynomials.

    Notes:

      Trailing zero terms are stripped away upon instantiation.

      Not None but, rather, p = Polynomial([]) is the 0 polynomial; for which p.degree()
      returns None.

      A non-zero constant polynomial has degree 0.

      Be careful: comparison is just comparison of tuples which may or may not be what
      the user wants on the level of polynomials; and, similarly, regarding max, min,
      sort, reverse, etc

      All polynomial operations (addition, muliplication, etc.) return instances of the
      same type as self; that is, a Polynomial or a possibly a subclass of Polynomial if
      in fact this class has been extended.

    Todo:

      1. Consider implementing improved typing in which the coefficients are subclasses
         of Number.
    """

    """Implementation notes:

      Some methods use local variables (which is more efficient than accessing instance
      variables) whenever, in general, doing so improves performance.

      class invariant:

      If poly is instantiated using poly = Polynomial([a_0,a_1,...,a_n]), then

        1. self._degree is the mathematical degree of poly:
              - if poly is not the 0 polynomial, then its degree is the largest
                m such that a_m is not zero.
              - the degree of 0 is None.

        2. len(self) = self._degree + 1

        3. self._coeffs is a tuple consisting of the coefficients of self.
    """

    def __init__(self, coeffs: Sequence[Ring] = ()):
        """Create a Polynomial.

        Polynomial([a_0,a_1,...,a_n]) or Polynomial((a_0,a_1,...,a_n)) constructs
        the polynomial a_0 + a_1 x + ... + a_n x^n.  In either case, the coeffici-
        ents will be an indexed tuple.

        Interactively, e.g.,

        >>> p = Polynomial([-1, 2, -3, 0, .4, 0])
        >>> print(p)
        -1 + 2x - 3x^2 + 0.4x^4
        >>> print(p.degree())
        4
        >>> p
        Polynomial((-1, 2, -3, 0, 0.4))
        >>> p[2]
        -3
        >>> print(p[1:])
        (2, -3, 0, 0.4)
        >>> print(Polynomial((0,) + p[1:]))
        2x - 3x^2 + 0.4x^4
        >>> print(Polynomial([0] + list(p[1:])))
        2x - 3x^2 + 0.4x^4

        >>> from fractions import Fraction
        >>> q = Polynomial((Fraction(1), Fraction(1, 3), Fraction('-2/5')))
        >>> print(q)
        1 + 1/3x - 2/5x^2
        >>> q
        Polynomial((Fraction(1, 1), Fraction(1, 3), Fraction(-2, 5)))
        >>> q = Polynomial((Fraction(1), Fraction('1/3'), Fraction(-0.4).limit_denominator()))
        >>> print(q)
        1 + 1/3x - 2/5x^2
        >>> q = Polynomial((Fraction(1), Fraction('1/3'), Fraction('-0.4')))
        >>> print(q)
        1 + 1/3x - 2/5x^2
        >>> from decimal import Decimal, getcontext
        >>> getcontext().prec = 5
        >>> q = Polynomial((Decimal(1), Decimal(1)/Decimal(3), Decimal('-0.4')))
        >>> print(q)
        1 + 0.33333x - 0.4x^2

        """
        """pre: coeffs is a sequence of coefficients for the polynomial being instantiated
           post: self._degree has been set to the mathematical degree of the polynomial
                 self._coeffs is a list of length self._degree the Polynomials coefficients
        """
        index = len(coeffs) - 1  # index is the index of the last nonzero coefficient
        while index > -1 and coeffs[index] == 0:
            coeffs = coeffs[:index]
            index -= 1
        if index == -1:
            self._coeffs: Sequence[Ring] = (0,)
            self._degree = None
        else:
            self._coeffs = tuple(coeffs)
            self._degree = index

    def __add__(self, other: Union[Number, Polynomial]):
        """Return the sum of two Polynomials.

        Coerces constants into constant Polynomials.

        >>> p1 = Polynomial([1, 2, 3])
        >>> print(p1+3)
        4 + 2x + 3x^2

        >>> p1 = Polynomial([1, 2, 3])
        >>> print(3+p1)
        4 + 2x + 3x^2

        >>> p2 = Polynomial([1/2, 2]);
        >>> print(p1+3+p2)
        4.5 + 4x + 3x^2

        """
        """pre: self and other's type is Polynomial or Number.
           post: return the sum of self and other as polynomials.
        """
        if isinstance(other, Ring):
            return self._add(self.__class__((other,)))
        elif isinstance(other, Polynomial):
            return self._add(other)
        else:
            return NotImplemented

    def __radd__(self, other: Union[Ring, Polynomial]):
        """Reverse add."""
        if isinstance(other, Ring):
            return self._add(self.__class__((other,)))
        elif isinstance(other, Polynomial):
            return self._add(other)
        else:
            return NotImplemented

    def _add(self, other: Polynomial):
        """Addition helper."""
        if self._degree is None:
            return other
        if other._degree is None:
            return self
        mindeg = min([self._degree, other._degree])
        if mindeg == self._degree:
            shorter = tuple(self._coeffs)
            longer = tuple(other._coeffs)
        else:
            shorter = tuple(other._coeffs)
            longer = tuple(self._coeffs)
        return self.__class__(
            tuple([shorter[i] + longer[i] for i in range(mindeg + 1)])
            + longer[mindeg + 1 :]
        )

    def __neg__(self):
        """Return the negative of a Polynomial.

        >>> p = Polynomial([0, 2, 3])
        >>> print(-p)
        -2x - 3x^2

        """
        return self.__class__([-x for x in self])

    def __sub__(self, other):
        """Return the difference of two Polynomials, interpreting scalars as Polynomials.

        Coerces constants into constant Polynomials.

        >>> p = Polynomial([0, 2, 3])
        >>> q = Polynomial([1, 2, 3])
        >>> print(p - q)
        -1

        >>> print(p - 1)
        -1 + 2x + 3x^2

        >>> p - 1
        Polynomial((-1, 2, 3))

        """
        """pre: self and other's type is Polynomial or Number.
           post: return the difference of self and other as polynomials."""
        return self + other.__neg__()

    def __mul__(self, other):
        """Return the product of two Polynomials, interpreting scalars as Polynomials.

        Coerces constants into constant Polynomials.

        >>> from fractions import Fraction
        >>> print(Polynomial([1, .2, Fraction(2, 3)]) * 2)
        2 + 0.4x + 4/3x^2

        >>> print(Polynomial([1, 2, Fraction(2, 3)]) * 2)
        2 + 4x + 4/3x^2

        >>> print(2 * Polynomial([1, 2, Fraction(2, 3)]))
        2 + 4x + 4/3x^2

        >>> print(Polynomial([1, Fraction(1, 2)]) * Polynomial([Fraction(1, 3), Fraction(1, 5)]))
        1/3 + 11/30x + 1/10x^2

        >>> print(Polynomial([1])*Polynomial([1,2]))
        1 + 2x

        """
        """pre: self and other's type is Polynomial or Number.
           post: return the sum of self and other as polynomials."""
        if isinstance(other, Number):
            return self._mul(self.__class__([other]))
        elif isinstance(other, Polynomial):
            return self._mul(other)
        else:
            return NotImplemented

    def __rmul__(self: Polynomial, other: Union[Ring, Polynomial]) -> Polynomial:
        """Reverse multiply."""
        if isinstance(other, Ring):
            return self._mul(self.__class__([other]))
        elif isinstance(other, Polynomial):
            return self._mul(other)
        else:
            return NotImplemented

    def _mul(self: Polynomial, other: Polynomial) -> Polynomial:
        """Multiplication helper."""
        if self._degree is None or other._degree is None:
            return self.__class__([0])
        product = []

        # See chapter 17, section 17.2, the section on vector convolutions in the
        # text Algorithms and Theory of Computation Handbook (1999) for the starting
        # point for deriving the algorithm below.
        lowerdeg = min([self._degree, other._degree])
        if self._degree == lowerdeg:
            shorter = self._coeffs
            longer = other._coeffs
            higherdeg = other._degree
        else:
            shorter = other._coeffs
            longer = self._coeffs
            higherdeg = self._degree
        for i in range(higherdeg + 1):
            summa = 0 # type: Ring
            if i <= lowerdeg:
                for j in range(i + 1):
                    summa += shorter[j] * longer[i - j]
                product.append(summa)
            else:
                for j in range(lowerdeg + 1):
                    summa += shorter[j] * longer[i - j]
                product.append(summa)
        for i in range(lowerdeg):
            summa_ = 0 # type: Ring
            for j in range(i + 1, lowerdeg + 1):
                summa_ += shorter[j] * longer[higherdeg + 1 + i - j]
            product.append(summa_)
        return self.__class__(product)

        ##Similar to the above algorithm but uses zero padding. Slightly slower,
        ##in general.
        # n = max([self._degree, other._degree])
        # if other._degree == n:
        #    lst1 = list(self) + (n - self._degree) * [0]
        #    lst2 = list(other)
        # else:
        #    lst1 = list(other) + (n - other._degree) * [0]
        #    lst2 = list(self)        # lst1 and lst2 both now have length n.
        # product = []
        # for i in range(n+1):         # compute the first n+1 coefficients of
        #    summ = 0                 # of the product.
        #    for j in range(i+1):
        #        summ += lst1[j] * lst2[i-j]
        #    product.append(summ)
        # for i in range(n):         # compute the rest of the coefficients of
        #    summ = 0                 # of the product of an n degree poly and
        #    for j in range(i+1,n+1): # a zero-padded n degree poly.
        #        summ += lst1[j] * lst2[n+1+i-j]
        #    product.append(summ)
        # return self.__class__(product)

    def fftmult(self: Polynomial, other: Polynomial) -> Polynomial:
        """Return the product of two polynomials, computed using (numpy's) FFT.

        >>> p = Polynomial([1,2,3]); q = Polynomial([1,2]);
        >>> print(p.fftmult(q))
        1.0 + 4.0x + 7.0x^2 + 6.0x^3

        Notes: O(nlog(n)). Doesn't work without modification with Fraction coefficients.
        """
        import numpy as np

        if self._degree is None or other._degree is None:
            return self.__class__([0])
        m = self._degree
        n = other._degree

        zero: List[Ring] = [0]
        p1 = np.array(list(self._coeffs) + n * zero)
        p2 = np.array(list(other._coeffs) + m * zero)
        return self.__class__((np.fft.ifft(np.fft.fft(p1) * np.fft.fft(p2))).real)

    def degree(self: Polynomial) -> Optional[int]:
        """Return the degree of a Polynomial.

        >>> p = Polynomial([1, 2, 0, 0])
        >>> print(p,"has degree", p.degree())
        1 + 2x has degree 1

        >>> p = Polynomial([17])
        >>> print(p.degree())
        0

        >>> p = Polynomial([])
        >>> p
        Polynomial((0,))
        >>> print(p.degree())
        None

        """
        return self._degree

    def __pow__(self: Polynomial, n: int) -> Polynomial:
        """Return the non-negative integral power of a Polynomial, computed recursively.

        >>> p = Polynomial([2])
        >>> print(p ** 6)
        64

        >>> print(Polynomial([3, 2])**2)
        9 + 12x + 4x^2

        >>> from fractions import Fraction
        >>> print(Polynomial([1, Fraction(1,2)])**2)
        1 + x + 1/4x^2

        >>> from fractions import Fraction
        >>> print(Polynomial([1, Fraction(1,2)])**3)
        1 + 3/2x + 3/4x^2 + 1/8x^3

        """

        def recpow(p, n):
            if n == 0:
                return self.__class__([1])
            else:
                factor = recpow(p, n // 2)
                # print("in recpow",factor)
                if n % 2 == 0:
                    return factor * factor
                else:
                    return factor * factor * p

        return recpow(self, n)

    def fftpow(self: Polynomial, n: int) -> Polynomial:
        """Return the non-negative integral power of a Polynomial.

        Computed recursively using fftmult.

        Note: This shouldn't be used, without modification, for high powers
        and may exhibit substantial errors even for low powers.
        >>> print(Polynomial([1, 2]).fftpow(3))
        1.0 + 6.0x + 12.0x^2 + 8.0x^3

        """

        def fftrecpow(p, n):
            if n == 0:
                return Polynomial([1])
            else:
                factor = fftrecpow(p, n // 2)
                # print("in recpow",factor)
                if n % 2 == 0:
                    return factor.fftmult(factor)
                else:
                    return factor.fftmult(factor.fftmult(p))

        return fftrecpow(self, n)

    def of(self: Polynomial, x: Ring) -> Ring:
        """Return the result of evaluating a Polynomial on a number.

        >>> p = Polynomial([1, 17, -1])
        >>> print(p.of(1))
        17
        >>> f = p.of
        >>> f(-1)
        -17
        >>> f(3)
        43

        """
        result: Ring = 0
        if self._degree is None:
            return result
        else:
            for i in range(self._degree + 1):
                result += self[i] * x ** i
            return result

    def __str__(self: Polynomial) -> str:
        """String coercion.

        >>> print(Polynomial([0,-1,2,0,-1,1,-2.3]))
        -x + 2x^2 - x^4 + x^5 - 2.3x^6

        """
        s = ""
        if self._degree == 0:
            return str(self[0])
        if self._degree is None:
            return "0"
        for i in range(0, self._degree + 1):
            if self[i] > 0:  # add coefficient
                if i != 0 and self[i] == 1 and s != "":
                    s += " + "
                elif i != 0 and s != "":
                    s += " + " + str(self[i])
                else:
                    s += str(self[i])
            elif self[i] < 0:
                if self[i] == -1 and s == "":
                    if i == 0:
                        s += "-1"
                    else:
                        s += "-"
                elif self[i] == -1:
                    s += " - "
                elif s == "":
                    s += "-" + str(-self[i])
                else:
                    s += " - " + str(-self[i])
            if i > 1 and self[i] != 0:  # add x^n
                s += "x^" + str(i)
            elif i == 1 and self[i] != 0:  # add x
                s += "x"
        return s

    def __len__(self: Polynomial) -> int:
        """Return number of coefficients of a Polynomial, which is its degree + 1."""

        if self._degree is None:
            return 1
        else:
            return self._degree + 1

    def __eq__(self: Polynomial, other: object) -> bool:
        """Return true if the two Polynomials have the same degree and equal coefficients.

        Coerces constants into constant Polynomials.

        >>> Polynomial([1]) == 1
        True

        >>> from fractions import Fraction
        >>> p1 = Polynomial([ 1, 2, 3, 0, 0]);
        >>> p2 = Polynomial([Fraction(1, 2), Fraction(1), Fraction(3, 2)])
        >>> print("Is",p1,"== 2 * (", p2,") ?", p1 == 2 * p2)
        Is 1 + 2x + 3x^2 == 2 * ( 1/2 + x + 3/2x^2 ) ? True

        """
        if isinstance(other, Polynomial):
            return self._coeffs == other._coeffs
        elif isinstance(other, Ring):
            return self._coeffs == Polynomial([other])._coeffs
        else:
            return NotImplemented

    def __repr__(self: Polynomial) -> str:
        """Return repr string.

        >>> Polynomial([])
        Polynomial((0,))

        >>> Polynomial([1, 2, 3, 0])
        Polynomial((1, 2, 3))

        """
        return "Polynomial(%s)" % str(self._coeffs)

    def __getitem__(self: Polynomial, n: int) -> Ring:
        """Built-in Indexing.

        Return an element tuple self._coeffs.

        >>> p = Polynomial((1, 2, 3, 4, 0))
        >>> print(p[2:])
        (3, 4)

        >>> it = iter(p)
        >>> next(it)
        1
        >>> next(it)
        2

        >>> print(p[-1])
        4

        """
        return self._coeffs[n]

    def __copy__(self: Polynomial) -> Polynomial:
        new_instance = type(self)(self._coeffs)
        new_instance.__dict__.update(self.__dict__)
        return new_instance

    def __deep__copy(self: Polynomial, memodict: dict = {}):
        new_instance = type(self)(copy.deepcopy(self._coeffs, memodict))
        new_instance.__dict__.update(self.__dict__)
        return new_instance


class FPolynomial(Polynomial):

    """Extends the Polynomial class for use with polynomials over a field.

    This subclass adds two methods, __divmod__ and __mod__, to the Polynomial class.
    each of which return polynomial(s) of type FPolynomial.

    """

    def divmod(self, other):

        """Return the quotient and remainder when dividing self by other.

        >>> from fractions import Fraction
        >>> numerator = FPolynomial((1,0,0,0,Fraction(-1)))
        >>> divisor = FPolynomial((1, Fraction(-2)))
        >>> print(numerator); print(divisor)
        1 - x^4
        1 - 2x
        >>> q, r = numerator.divmod(divisor)
        >>> print(q); print(r)
        1/16 + 1/8x + 1/4x^2 + 1/2x^3
        15/16
        """
        if isinstance(other, Number):
            other = FPolynomial([other])
        elif not isinstance(other, FPolynomial):
            return NotImplemented

        if self.degree() is None:
            return self.__class__([0]), self.__class__([0])

        assert other.degree is not None, "Cannot divide by zero."

        num = FPolynomial(self._coeffs)
        quo = FPolynomial([])

        if num.degree() >= other.degree():
            while num.degree() is not None and num.degree() >= other.degree():
                monomial = FPolynomial(
                    (num.degree() - other.degree()) * [0] + [num[-1] / other[-1]]
                )
                num = num - monomial * other
                quo = quo + monomial
            return self.__class__(quo._coeffs), self.__class__(num._coeffs)
        else:
            return self.__class__([0]), self

    def __mod__(self, other):

        """Return the remainder when dividing self by other.

        >>> numerator = FPolynomial((1,0,0,0,-1))
        >>> divisor = FPolynomial((1,-1))
        >>> print(numerator); print(divisor)
        1 - x^4
        1 - x
        >>> print(numerator % divisor)
        0
        """
        if isinstance(other, Number):
            other = FPolynomial([other])
        elif not isinstance(other, FPolynomial):
            return NotImplemented

        if self.degree() is None:
            return FPolynomial([0])

        assert other.degree is not None, "Cannot divide by zero."

        num = FPolynomial(self._coeffs)

        if num.degree() >= other.degree():
            while num.degree() is not None and num.degree() >= other.degree():
                monomial = FPolynomial(
                    (num.degree() - other.degree()) * [0] + [num[-1] / other[-1]]
                )
                num = num - monomial * other
            return self.__class__(num._coeffs)
        else:
            return self

    def __repr__(self):
        """Return repr string.

        >>> FPolynomial([])
        FPolynomial((0,))

        >>> FPolynomial([1, 2, 3, 0])
        FPolynomial((1, 2, 3))
        """
        return "FPolynomial(%s)" % str(self._coeffs)


if __name__ == "__main__":

    import doctest

    doctest.testmod()
