"""Provide polynomial classes.

 The base class here is a polynomial class implemented using composition with
 tuple objects.  The underlying concrete representation of a polynomial is a
 tuple of its coefficients of length one less than the polynomial's degree.
 (The zero polynomial has degree None.)

 To see several examples at your commandline:

     pydoc polylib.Polynomial.__init__

 Also, c.f.

     pydoc polylib.FPolynomial

 Notes on typing:

 This is type-hinted throughout (which does not provide mammoth utility). To
 demonstrate what one can in fact do, if you create the following program and
 call it test_typing.py:

     from polylib import Polynomial
     p = Polynomial(['adsf','ghj'])

 and then type check it with:

     mypy test_typing.py

 you'll get an error to the effect:

     Value of type variable "Ring" of "Polynomial" cannot be "str"

 This happens because Python's str class (which does implement addition) does
 not implement, for example, multiplication (of strings) and so can't be used
 for a ring.

 However, note that, even though integers are not a division ring, one does
 not, unfortunately, get a type error with

     p = FPolynomial([1, 2, 3]) # type is polynomial.FPolynomial[builtin.int*]

 this is due to the fact that Python's int implements (true) division (return-
 ing a float).  (Currently, mypy's Protocol types just check for the existence
 of a __truediv__ method not its type signature.)

 So, checks involving typing are limited.

 Also, all rings are currently assumed to be commutative but may work correcty
 with noncommutative rings - it's just that that hasn't been thoroughly checked.

 TODO:

 - When we can assume Python10, no need to import annotations from __future__
   (which is currently used for the type hints for other).
 - If/when Python's Number numeric hierarchy works with type hints, one might
   be able to use nominal typing instead of structural typing.
   Or, or might be better to keep structural typing in any case (so we can use
   gmpy2 or numpy ints or any numbers that cannot be made to subclass Number).
   UPDATE: gmpy2 and numpy types may subclass Number.
 - Later versions of mypy import Protocol from typing, not typing_extensions.
   Same for runtime which may be called runtime_checkable in the future.
 - Verify that this works with coefficients in a non-cummutative ring and/or
   tweak the dunder methods so that it does work.
"""

from __future__ import annotations
import copy
from typing import Sequence, Union, Tuple, List, Optional, cast, TypeVar, Generic
import sys
if sys.version_info >= (3, 8):
    from typing import Protocol
else:
    from typing_extensions import Protocol
try:
    from typing import runtime_checkable as runtime
except:
    from typing_extensions import runtime

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


@runtime
class Ring_(Protocol):

    """The minimal operations that we require for a ring.

    This in only used for structural typing purposes.
    """

    def __add__(self, other):
        ...

    def __neg__(self):
        ...

    def __sub__(self, other):
        ...

    def __mul__(self, other):
        ...

    def __pow__(self, n):
        ...


@runtime
class OrderedRing_(Ring_, Protocol):

    """A totally ordered ring."""

    def __gt__(self, other) -> bool:
        ...

    def __lt__(self, other) -> bool:
        ...


@runtime
class DivisionRing_(Ring_, Protocol):

    """A ring invertible non-zero elements (possibly noncommutative)."""

    def __truediv__(self, other):
        ...


Ring = TypeVar("Ring", bound="Ring_")

# Below won't work because currently a TypeVar can't take an argument. If this worked we
# might be able to put P[R] instead Polynomial[R] in the dunder methods in Polynomial.
# Then, for example, we wouldn't have to cast to FPolynomial in the FPolynomial class.
# The issue is that Python typing can't figure out subclasses of Polynomial.
# P = TypeVar('P', bound='Polynomial')


class Polynomial(Generic[Ring]):

    """Implements polynomials over a ring.

    Addition, subtraction, multiplication, and evaluation of polynomials is
    implemented via dunder methods.

    Notes:

      The coefficients can be in any implementation of a ring, possibly non-
      commutative and without 1, though

         - said implementation of coefficients must coerce (right) addition
           of an element by the int 0.
         - if __truediv__ is implemented then we assume (only in __str__)
           that the ring has a unit (and that right multiplication by 1 is
           implemented).
         - TODO: define a non-commutative ring and verify that rmul, etc.
           is not in reality assuming commutativity (and fix, if necessary).

      Trailing zero terms of polynomials are stripped away upon instantiation.

      Polynomial([0]), where 0 is the zero from the coefficient ring, is the
      zero polynomial; for which p.degree() returns, and p._degree is, None.

      A non-zero constant polynomial has degree 0.

      The argument to Polynomial() should be a non-empty Sequence such as a
      list or a tuple.

      Be careful: __eq__ is just comparison of tuples which may not be what
      the user wants on the level of polynomials; and, similarly, regarding
      max, min, sort, reverse, etc.

      All polynomial operations (addition, multiplication, etc.) return in-
      stances of the same type as self; that is, a Polynomial or a possibly
      a subclass of Polynomial if in fact this class has been extended.
    """

    """Implementation notes:

      Some methods use local variables (which is more efficient than access-
      ing instance variables) whenever doing so improves performance.

      class invariant:

      If poly is instantiated using poly = Polynomial([a_0,a_1,...,a_n]),
      then

        1. self._degree is the mathematical degree of poly:
              - if poly is not the 0 polynomial, then its degree is the
                largest m such that a_m is not zero.
              - the degree of 0 is None.

        2. len(self) = self._degree + 1

        3. self._coeffs is a tuple consisting of the coefficients of self.
    """
    __slots__ = ('_coeffs', '_degree', 'x', 'spaces', 'increasing', 'x_unwrapped')

    def __init__(self: Polynomial[Ring], coeffs: Sequence[Ring], x: str = "x", spaces: bool = True, increasing: bool = True) -> None:
        """Create a Polynomial.

        Polynomial([a_0,a_1,...,a_n]) or Polynomial((a_0,a_1,...,a_n)) con-
        structs the polynomial a_0 + a_1 x + ... + a_n x^n.  In either case,
        the coefficients will be an indexed tuple.

        The argument to the parameter x controls the letter used for the in-
        determinant and whether the str representation wraps the polynomial
        in parentheses or in square brackets in which case, with a totally
        ordered coefficient ring, the value of spaces is automatically set
        to False (see examples).

        Args:

            coeffs (iterable): The coefficients for the polynomial.

            x (string): The string to use for the indeterminate.

            spaces (bool): True leads to removing spaces from the string
                representation of polynomials, but only when the coeffi-
                cient ring is totally ordered.

            increasing (bool): Whether the string representing has inc-
                creasing powers of the indeterminant.

        Examples:

            >>> p = Polynomial([-1, 2, -3, 0, 4, 0])
            >>> print(p)
            -1 + 2x - 3x^2 + 4x^4
            >>> print(p.degree())
            4
            >>> p
            Polynomial((-1, 2, -3, 0, 4))
            >>> p[2]
            -3
            >>> p[1:] # Note that this is not an instance of Polynomial
            (2, -3, 0, 4)

            >>> Polynomial(())
            Traceback (most recent call last):
            ValueError: coeffs cannot be empty

            >>> p = Polynomial([3]); p
            Polynomial((3,))
            >>> print(p.degree())
            0

            >>> p = Polynomial([0]); p
            Polynomial((0,))
            >>> print(p.degree())
            None

            >>> from fractions import Fraction as Frac
            >>> q = Polynomial((Frac(1), Frac(1, 3), Frac('-2/5')))
            >>> print(q)
            1 + 1/3x - 2/5x^2
            >>> q
            Polynomial((Fraction(1, 1), Fraction(1, 3), Fraction(-2, 5)))
            >>> q = Polynomial((1, Frac('1/3'), Frac(-0.4).limit_denominator()))
            >>> print(q)
            1 + 1/3x - 2/5x^2
            >>> q = Polynomial((Frac(1), Frac('1/3'), Frac('-0.4')))
            >>> print(q)
            1 + 1/3x - 2/5x^2

            >>> from decimal import Decimal, getcontext
            >>> getcontext().prec = 5
            >>> q=Polynomial((1,Decimal(1)/Decimal(3),Decimal('-0.4')))
            >>> print(q)
            1 + 0.33333x - 0.4x^2

            >>> t = Polynomial([0, 1], x='t')
            >>> print((3*t**2-5))
            -5 + 3t^2

            >>> t = Polynomial([0, 1], x='t', increasing = False)
            >>> print((3*t**2-5))
            3t^2 - 5

            >>> print(Polynomial([4*t**0-5*t, 3*t**3+t**4]))
            -5t + 4 + t^4 + 3t^3x
            >>> # str representation on the last line is ambiguous
            >>> t = Polynomial([0, 1], x='(t)')   # wrap the coefficient polys
            >>> print(Polynomial([4*t**0-5*t, 3*t**3+t**4]))  # use t**0 for 1
            (4-5t) + (3t^3+t^4)x

            >>> t = Polynomial([0, 1], x='[t]', increasing = False)
            >>> p=Polynomial([complex(4)*t**0-complex(5,4)*t, complex(3)*t**3])
            >>> print(p)
            [(-5-4j)t + (4+0j)] + [(3+0j)t^3]x

            >>> t = Polynomial([0, 1], x='[t]')
            >>> p=Polynomial([complex(4)*t**0-complex(5,4)*t,complex(3)*t**3],'X')
            >>> print(p)
            [(4+0j) + (-5-4j)t] + [(3+0j)t^3]X

            >>> x = Polynomial([0, 1])
            >>> t = Polynomial([0, 1], '(t)')
            >>> print(t * x ** 0 + (t**2) * x)
            (t) + (t^2)x
            >>> print(x + t)
            (t) + x
            >>> print(t + x)
            (x + t)
            >>> x + t == t + x
            True

            >>> p = (t**0 + 2*t)*x**0 + (3*t+4*t**2)*x
            >>> p
            Polynomial((Polynomial((1, 2)), Polynomial((0, 3, 4))))
            >>> print(p)
            (1+2t) + (3t+4t^2)x
        """
        if len(coeffs) == 0:
            raise ValueError("coeffs cannot be empty")

        self.x = x
        self.increasing = increasing

        if (x[0] == "(" and x[-1] == ")") or (x[0] == "[" and x[-1] == "]"):
            self.x_unwrapped = x[1:-1]
            self.spaces = False
        else:
            self.x_unwrapped = x
            self.spaces = spaces

        self._degree: Optional[int]

        index = len(coeffs) - 1  # will be the index of the last nonzero coeff
        while index > 0 and coeffs[index] == 0:  # remove zeros
            coeffs = coeffs[:index]
            #coeffs = coeffs[:-1]
            index -= 1

        # slower
        #index = len(coeffs) - 1  # will be the index of the last nonzero coeff
        #while index > 0 and coeffs[index] == 0:  # remove zeros
        #    index -= 1
        #coeffs = coeffs[:index+1]

        # about same as first way
        #while len(coeffs) > 1 and coeffs[-1] == 0:
        #      coeffs = coeffs[:-1]
        #index = len(coeffs) - 1

        if index > 0:
            self._degree = index
        else:
            self._degree = None if coeffs[0] == coeffs[0] - coeffs[0] else 0
        #coeffs = [
        #    cast(Ring, 0) * coeffs[-1] if coeff == 0 else coeff for coeff in coeffs
        #]
        self._coeffs: Sequence[Ring] = tuple(coeffs)

    def __add__(
        self: Polynomial[Ring], other: Union[Ring, Polynomial[Ring]]
    ) -> Polynomial[Ring]:
        """Return the sum of two Polynomials.

        Coerces constants into constant Polynomials.

        Examples:

            >>> p1 = Polynomial([1, 2, 3])
            >>> print(p1+3)
            4 + 2x + 3x^2

            >>> print(3+p1)
            4 + 2x + 3x^2

            >>> p2 = Polynomial([1/2, 2]);
            >>> print(p1+3+p2)
            4.5 + 4x + 3x^2
        """
        if isinstance(other, Polynomial) and self.x_unwrapped == other.x_unwrapped:
            selfdeg = self._degree; otherdeg = other._degree
            selfcos = self._coeffs; othercos = other._coeffs
            if selfdeg is None:
                return other
            if otherdeg is None:
                return self
            mindeg = min([selfdeg, otherdeg])
            if mindeg == selfdeg:
                return self.__class__(
                    tuple(selfcos[i] + othercos[i] for i in range(mindeg + 1)) + othercos[mindeg + 1 :],
                    self.x,
                    self.spaces,
                    self.increasing
                )
            else:
                return self.__class__(
                    tuple(selfcos[i] + othercos[i] for i in range(mindeg + 1)) + selfcos[mindeg + 1 :],
                    self.x,
                    self.spaces,
                    self.increasing
                )
        return self.__class__((self[0] + other,)+self[1:], self.x, self.spaces, self.increasing)

    def __radd__(self: Polynomial[Ring], other: Ring) -> Polynomial[Ring]:
        """Reverse add."""
        return self.__class__((self[0] + other,)+self[1:], self.x, self.spaces, self.increasing)

    def __neg__(self: Polynomial[Ring]) -> Polynomial[Ring]:
        """Return the negative of a Polynomial.

        Examples:

            >>> p = Polynomial([0, 2, 3])
            >>> print(-p)
            -2x - 3x^2
        """
        return self.__class__([-co for co in self._coeffs], self.x, self.spaces, self.increasing)

    def __sub__(self: Polynomial[Ring], other: Union[Ring, Polynomial[Ring]]) -> Polynomial[Ring]:
        """Return the difference of two Polynomials.

        Coerces constants into constant Polynomials.

        Examples:

            >>> p = Polynomial([0, 2, 3])
            >>> q = Polynomial([1, 2, 3])
            >>> print(p - q)
            -1

            >>> print(p - 1)
            -1 + 2x + 3x^2

            >>> p - 1
            Polynomial((-1, 2, 3))
        """
        return self + other.__neg__()

    def __rsub__(self: Polynomial[Ring], other: Ring) -> Polynomial[Ring]:
        """Reverse subtract."""
        #if isinstance(other, Ring_):
        #return (- self).__add__(self.__class__((other,), self.x, self.spaces, self.increasing))
        return self.__class__((-self[0] + other,)+tuple(-co for co in self._coeffs[1:]), self.x, self.spaces, self.increasing)
        #return NotImplemented

    def __mul__(
        self: Polynomial[Ring], other: Union[Ring, Polynomial[Ring]]
    ) -> Polynomial[Ring]:
        """Return the product of two Polynomials.

        Coerces constants into constant Polynomials.

        Examples:

            >>> from fractions import Fraction
            >>> print(Polynomial([1, .2, Fraction(2, 3)]) * 2)
            2 + 0.4x + 4/3x^2

            >>> print(Polynomial([1, 2, Fraction(2, 3)]) * 2)
            2 + 4x + 4/3x^2

            >>> print(2 * Polynomial([1, 2, Fraction(2, 3)]))
            2 + 4x + 4/3x^2

            >>> p1 = Polynomial([1, Fraction(1, 2)])
            >>> p2 = Polynomial([Fraction(1, 3), Fraction(1, 5)])
            >>> print(p1 * p2)
            1/3 + 11/30x + 1/10x^2
        """
        if isinstance(other, Polynomial):
            if self.x_unwrapped == other.x_unwrapped:
                selfdeg = self._degree; otherdeg = other._degree
                selfcos = self._coeffs; othercos = other._coeffs
                if selfdeg is None or otherdeg is None:
                    return self.__class__([cast(Ring, 0)], self.x, self.spaces, self.increasing)
                product: List[Ring] = []
                # See chapter 17, section 17.2, the section on vector convolutions in the
                # text Algorithms and Theory of Computation Handbook (1999) for the starting
                # point for deriving the algorithm below.
                lowerdeg = min([selfdeg, otherdeg])
                if selfdeg == lowerdeg:
                    shorter = selfcos
                    longer = othercos
                    higherdeg = otherdeg
                else:
                    shorter = othercos
                    longer = selfcos
                    higherdeg = selfdeg
                for i in range(higherdeg + 1):
                    summa = cast(Ring, 0)
                    if i <= lowerdeg:
                        for j in range(i + 1):
                            summa = shorter[j] * longer[i - j] + summa
                        product.append(summa)
                    else:
                        for j in range(lowerdeg + 1):
                            summa = shorter[j] * longer[i - j] + summa
                        product.append(summa)
                for i in range(lowerdeg):
                    summa_ = cast(Ring, 0)
                    for j in range(i + 1, lowerdeg + 1):
                        summa_ = shorter[j] * longer[higherdeg + 1 + i - j] + summa_
                    product.append(summa_)
                return self.__class__(product, self.x, self.spaces, self.increasing)

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
                # return self.__class__(product, self.x, self.spaces, self.increasing)
            if other[-1] != other[-1]**0 and other[:-1] != 0:
                return NotImplemented
            return other.__class__(other._degree * (0,) + (self,), other.x)
            #return self.__class__(tuple(self * coef for coef in other._coeffs), other.x, other.spaces, other.increasing)
        #if isinstance(other, Ring_):
        #return self._mul(self.__class__((other,), self.x, self.spaces, self.increasing))
        return self.__class__(tuple(coef * other for coef in self._coeffs), self.x, self.spaces, self.increasing)
        #return NotImplemented

    def __rmul__(self: Polynomial[Ring], other: Ring) -> Polynomial[Ring]:
        """Reverse multiply."""

        #if isinstance(other, Ring_):
            #return self.__class__((other,), self.x, self.spaces, self.increasing)._mul(self)
        #return self._mul(self.__class__((other,), self.x, self.spaces, self.increasing))
        return self.__class__(tuple(coef * other for coef in self._coeffs), self.x, self.spaces, self.increasing)
        #return NotImplemented

    def fftmult(self: Polynomial[Ring], other: Polynomial[Ring]) -> Polynomial[Ring]:
        """Return the product of two polynomials, computed using (numpy's) FFT.

        >>> p = Polynomial([1,2,3]); q = Polynomial([1,2]);
        >>> print(p.fftmult(q))
        1.0 + 4.0x + 7.0x^2 + 6.0x^3

        Notes: O(nlog(n)). Doesn't work without modification with Fraction
        coefficients.
        """
        import numpy as np

        if self._degree is None or other._degree is None:
            return self.__class__([cast(Ring, 0)], self.x, self.spaces, self.increasing)
        m = self._degree
        n = other._degree

        zero: List[Ring] = [cast(Ring, 0)]
        p1 = np.array(list(self._coeffs) + n * zero)
        p2 = np.array(list(other._coeffs) + m * zero)
        return self.__class__(
            (np.fft.ifft(np.fft.fft(p1) * np.fft.fft(p2))).real, self.x, self.spaces, self.increasing
        )

    def degree(self: Polynomial[Ring]) -> Optional[int]:
        """Return the degree of a Polynomial.

        Examples:

            >>> p = Polynomial([1, 2, 0, 0])
            >>> print(p, "has degree", p.degree())
            1 + 2x has degree 1

            >>> p = Polynomial([17])
            >>> print(p.degree())
            0

            >>> p = Polynomial([0])
            >>> p
            Polynomial((0,))
            >>> print(p.degree())
            None

            >>> p = Polynomial([complex(0)])
            >>> p
            Polynomial((0j,))
            >>> print(p.degree())
            None
        """
        return self._degree

    def __pow__(self: Polynomial[Ring], n: int) -> Polynomial[Ring]:
        """Return the non-negative integral power of a Polynomial.

        Examples:

            >>> p = Polynomial([2])
            >>> print(p ** 6)
            64

            >>> print(Polynomial([3, 2])**2)
            9 + 12x + 4x^2

            >>> Polynomial([3, 2])**0
            Polynomial((1,))

            >>> from fractions import Fraction
            >>> print(Polynomial([1, Fraction(1,2)])**2)
            1 + x + 1/4x^2

            >>> from fractions import Fraction
            >>> print(Polynomial([1, Fraction(1,2)])**3)
            1 + 3/2x + 3/4x^2 + 1/8x^3

            >>> p = Polynomial([0])
            >>> p**0
            Polynomial((1,))

            >>> p = Polynomial([1, 2, 7])
            >>> p**0
            Polynomial((1,))
        """
        if n < 0:
            raise ValueError(f"n must be nonnegative, not {n}")

        if self.degree() is None:
            if n == 0:
                return self.__class__([cast(Ring, 0) * self[-1] + cast(Ring, 1)], self.x, self.spaces, self.increasing)
            else:
                return self

        # NOTE: DO YOU WANT THIS OR JUST RECURSIVELY CALL POW
        def recpow(p: Polynomial[Ring], n: int) -> Polynomial[Ring]:
            if n == 0:
                return self.__class__(
                    [self[-1] * cast(Ring, 0) + cast(Ring, 1)], self.x, self.spaces, self.increasing
                )
            else:
                factor = recpow(p, n // 2)
                #print("in recpow", factor)
                if n % 2 == 0:
                    return factor * factor
                else:
                    return factor * factor * p

        return recpow(self, n)

    def fftpow(self: Polynomial[Ring], n: int) -> Polynomial[Ring]:
        """Return the non-negative integral power of a Polynomial.

        Computed recursively using fftmult.

        Note: This shouldn't be used, without modification, for high powers
        and may exhibit substantial errors even for low powers.

        Example:

            >>> print(Polynomial([1, 2]).fftpow(3))
            1.0 + 6.0x + 12.0x^2 + 8.0x^3
        """

        def fftrecpow(p, n):
            if n == 0:
                return Polynomial([1], self.x, self.spaces, self.increasing)
            else:
                factor = fftrecpow(p, n // 2)
                # print("in recpow",factor)
                if n % 2 == 0:
                    return factor.fftmult(factor)
                else:
                    return factor.fftmult(factor.fftmult(p))

        return fftrecpow(self, n)

    def __call__(self: Polynomial[Ring], x: Ring) -> Ring:
        """Return the result of evaluating a Polynomial on a number.

        Examples:

            >>> p = Polynomial([1, 17, -1])
            >>> print(p(1))
            17
            >>> p(-1)
            -17
            >>> p(3)
            43

            >>> Polynomial([0])(3)
            0

            >>> Polynomial((complex(0),))(3)
            0j

            We can also compose polynomials:
            >>> p = Polynomial([1, 2])  # 1 + 2x
            >>> print(p(p))
            3 + 4x
        """
        if self._degree is None or self.degree == 0 or x == 0 * self._coeffs[0]:
            return self._coeffs[0]
        result =  x ** self._degree * self._coeffs[-1]
        for i in range(self._degree):
            result += x ** i * self[i]
        return result

    def __str__(self: Polynomial[Ring], streamline: bool = True) -> str:
        """String coercion.

        Args:

          streamline (bool): True leads to not printing terms that have
              coefficient zero or one.

        Examples:

            >>> print(Polynomial([1, 2, 0, 0]))
            1 + 2x

            >>> print(Polynomial([0, 1]))
            x

            >>> p = Polynomial([1, 2, 0, -3])
            >>> print(p.__str__(streamline = False))
            1 + 2x + 0x^2 + -3x^3

            >>> p = Polynomial([1, 2, 0, -3], increasing = False)
            >>> print(p.__str__(streamline = False))
            -3x^3 + 0x^2 + 2x + 1

            >>> print(Polynomial([0, -1,2, 0, -1,1, -2.3]))
            -x + 2x^2 - x^4 + x^5 - 2.3x^6

            >>> p = Polynomial([complex(0), complex(-1,2), complex(0,-2.3)])
            >>> print(p)
            (-1+2j)x + -2.3jx^2

            >>> p = Polynomial([complex(1), complex(-1,2), complex(0,-2.3)])
            >>> print(p)
            1 + (-1+2j)x + -2.3jx^2

            >>> print(Polynomial([complex(1)]))
            (1+0j)

            >>> print(Polynomial([0, complex(1)]))
            x

            >>> print(Polynomial([complex(1), complex(-1,2), complex(1)]))
            1 + (-1+2j)x + x^2

            >>> p = Polynomial([ complex(-1,2), complex(-1,2), complex(1)])
            >>> print(p)
            (-1+2j) + (-1+2j)x + x^2

            Polynomials over Polynomials

            >>> t = Polynomial([0,1], x='t')
            >>> print(Polynomial([t**2+1,t-9]).__str__(streamline = False))
            1 + t^2 + -9 + tx

            >>> t = Polynomial([0,1], x='(t)')
            >>> print(Polynomial([t**2+1,t-9]).__str__(streamline = False))
            (1+t^2) + (-9+t)x

            >>> t = Polynomial([0,1], x='t')
            >>> print(Polynomial([t**2+1,t-9]))
            1 + t^2 + -9 + tx

            >>> t = Polynomial([0,1], x='(t)', spaces = False)
            >>> print(Polynomial([t**2+1, 9-t]))
            (1+t^2) + (9-t)x
            >>> print(Polynomial([t**2+1,t-9])**2)
            (1+2t^2+t^4) + (-18+2t-18t^2+2t^3)x + (81-18t+t^2)x^2
        """
        if self.x == self.x_unwrapped:
            lp = ""; rp = ""
            var = self.x
        else:
            lp = self.x[0]; rp = self.x[-1]
            var = self.x[1:-1]

        if self._degree == 0 or self._degree is None:
            return lp + str(self[0]) + rp
        s = ''
        def reverse(s):
            symbols = s.split()
            monomials = symbols[::-2]
            operators = symbols[-2::-2] + ['+']
            s = ' '.join(x + ' ' + y for x, y in zip(operators, monomials))
            s = s.replace('+ -','- ')
            s = s.replace('- -','+ ')
            if s[:2] == '+ ':
                s = s[2:]
            if s[:2] == '- ':
                s = '-' + s[2:]
            return s
        if streamline:
            elt: Ring_ = self._coeffs[-1]
            zero_: Ring_ = elt * 0
            if not hasattr(elt, "__truediv__"):
                streamline = False
            else:
                elt_: DivisionRing_ = elt + 1 if elt == zero_ else elt
                one: DivisionRing_ = elt_ / elt_
                #try:   # NOTE:  clean this up
                #    one: DivisionRing_ = elt_ / elt_
                #except:
                #    one: DivisionRing_ = elt_.__class__(1)
                try:  # successful if coeffs are OrderedRing
                    elt_ < elt_  # type: ignore
                    zero: OrderedRing_ = zero_ * elt
                    for i in range(0, self._degree + 1):
                        if self[i] > zero:  # add coefficient
                            if (
                                i != 0
                                and self[i] == one
                                and s != ""
                                and s != "("
                                and s != "["
                            ):
                                s += " + "
                            elif i != 0 and s != "" and s != "(" and s != "[":
                                s += " + " + str(self[i])
                            elif i == 0 or self[i] != one:
                                s += str(self[i])
                        elif self[i] < zero:
                            if self[i] == one.__neg__() and (
                                s == "" or s == "(" or s == "["
                            ):
                                if i == 0:
                                    s += "-1"
                                else:
                                    s += "-"
                            elif self[i] == one.__neg__():
                                s += " - "
                            elif s == "" or s == "(" or s == "[":
                                s += "-" + str(-self[i])
                            else:
                                s += " - " + str(-self[i])
                        if i > 1 and self[i] != zero:  # add x^n
                            s += var + "^" + str(i)
                        elif i == 1 and self[i] != zero:  # add x
                            s += var
                    if not self.increasing:
                        s = reverse(s)
                    if not self.spaces:
                        s = "".join(s.split())
                except:
                    zero__: DivisionRing_ = zero_ * elt
                    for i in range(0, self._degree + 1):
                        if i == 0 and self[0] != zero__:
                            if self[0] == one:
                                s += str(1) + " + "
                            else:
                                s += str(self[0]) + " + "
                        elif i > 1 and self[i] != zero__:  # add x^n
                            if self[i] == one:
                                s += var + "^" + str(i)
                            elif self[i] == -one:
                                s += "-" + var + "^" + str(i)
                            else:
                                s += str(self[i]) + var + "^" + str(i)
                            if i < self._degree:
                                s += " + "
                        elif i == 1 and self[i] != zero__:  # add x
                            if self[i] == one:
                                s += var
                            elif self[i] == -one:
                                s += "-" + var + "^" + str(i)
                            else:
                                s += str(self[i]) + var
                            if i < self._degree:
                                s += " + "
                    if not self.increasing:
                        s = reverse(s)
        if not streamline:
            for i in range(0, self._degree + 1):
                if i == 0:
                    s += str(self[0])
                    if i < self._degree:
                        s += " + "
                elif i > 1:  # add x^n
                    s += str(self[i]) + var + "^" + str(i)
                    if i < self._degree:
                        s += " + "
                elif i == 1:  # add x
                    s += str(self[i]) + var
                    if i < self._degree:
                        s += " + "
            if not self.increasing:
                s = reverse(s)
        return lp + s + rp

    def __len__(self: Polynomial[Ring]) -> int:
        """Return number of coefficients of a Polynomial, which is its degree + 1."""
        if self._degree is None:
            return 1
        else:
            return self._degree + 1

    def __eq__(self: Polynomial[Ring], other: object) -> bool:
        """Return true if the two Polynomials the same coefficients.

        Coerces constants into constant Polynomials.

        Examples:

            >>> Polynomial([1]) == 1
            True

            >>> from fractions import Fraction
            >>> p1 = Polynomial([ 1, 2, 3, 0, 0]);
            >>> p2 = Polynomial([Fraction(1, 2), Fraction(1), Fraction(3, 2)])
            >>> print("Is",p1,"== 2 * (", p2,") ?", p1 == 2 * p2)
            Is 1 + 2x + 3x^2 == 2 * ( 1/2 + x + 3/2x^2 ) ? True
        """
        if isinstance(other, Polynomial):
            return self._coeffs == other._coeffs and self.__class__ == other.__class__
        if isinstance(other, Ring_):
            return self._coeffs == self.__class__((other,))._coeffs
        return NotImplemented

    def __repr__(self: Polynomial[Ring]) -> str:
        """Return repr string.

        Examples:
            >>> Polynomial([0])
            Polynomial((0,))

            >>> Polynomial([complex(0)])
            Polynomial((0j,))

            >>> Polynomial([1, 2, 3, 0])
            Polynomial((1, 2, 3))
        """
        return "Polynomial(%s)" % repr(self._coeffs)

    def __getitem__(self: Polynomial[Ring], n: int) -> Ring:
        """Built-in Indexing.

        Return an element tuple self._coeffs.

        Examples:

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

    #def __copy__(self: Polynomial[Ring]) -> Polynomial[Ring]:
    #    new_instance = type(self)(self._coeffs)
    #    #new_instance.__dict__.update(self.__dict__)
    #    new_instance.__slots__ = self.__slots__
    #    return new_instance

    #def __deepcopy__(self: Polynomial[Ring], memodict: dict = {}) -> Polynomial[Ring]:
    #    new_instance = type(self)(copy.deepcopy(self._coeffs, memodict))
    #    new_instance.__dict__.update(self.__dict__)
    #    return new_instance

    def divmod(
        self: Polynomial[Ring], other: Polynomial[Ring]
    ) -> Tuple[Polynomial[Ring], ...]:
        """Return the quotient and remainder when dividing self by other.

        Examples:

            >>> numerator = Polynomial((5, 4, 3, 2, -1))
            >>> divisor = Polynomial((7, -1))
            >>> print(numerator); print(divisor)
            5 + 4x + 3x^2 + 2x^3 - x^4
            7 - x
            >>> q, r = numerator.divmod(divisor)
            >>> print(q); print(r)
            220 + 32x + 5x^2 + x^3
            -1535
            >>> numerator == q * divisor + r
            True
        """
        if not (isinstance(other, Polynomial) and other.x_unwrapped == self.x_unwrapped):
            #raise ValueError(f"{other.x_unwrapped} != {self.x_unwrapped}")
            return NotImplemented

        otherdeg = other._degree

        if otherdeg is None:
            raise ValueError("cannot divide by zero")

        # debug some stuff
        has_to_be_zero = (self[-1] - self[-1] * other[-1]) * (self[-1] + self[-1] * other[-1])
        if isinstance(has_to_be_zero, Polynomial):
            assert has_to_be_zero._degree is None, (
                f"potential infinite loop in long division in divmod neither "
                f"{self[-1]-self[-1]*other[-1]} nor {self[-1]+self[-1]*other[-1]} "
                f"are the zero poly, in Polynmial.divmod()"
            )
        else:
            assert has_to_be_zero == cast(Ring, 0) * self[-1], (
                f"potential infinite loop in long division in divmod neither "
                f"{self[-1]-self[-1]*other[-1]} nor {self[-1]+self[-1]*other[-1]} are "
                f"zero, in Polynmial.divmod()"
            )

        if self._degree is None:
            return (
                self.__class__([cast(Ring, 0)], self.x, self.spaces, self.increasing),
                self.__class__([cast(Ring, 0)], self.x, self.spaces, self.increasing),
            )

        num = copy.copy(self)
        quo = Polynomial([cast(Ring, 0)], self.x, self.spaces, self.increasing)

        numdeg = num._degree

        if cast(int, numdeg) >= cast(int, otherdeg):
            while numdeg is not None and cast(int, numdeg) >= cast( int, otherdeg):
                monomial = Polynomial(
                    (cast(int, numdeg) - cast(int, otherdeg))
                    * [cast(Ring, 0)]
                    + [num[-1] * other[-1]], self.x, self.spaces, self.increasing
                )
                num = cast(
                    Polynomial, num - monomial * other
                )  # shouldn't have to cast here?
                quo = cast(Polynomial, quo + monomial)
                numdeg = num._degree
            return (
                self.__class__(quo._coeffs, self.x, self.spaces, self.increasing),
                self.__class__(num._coeffs, self.x, self.spaces, self.increasing),
            )
        else:
            return self.__class__([cast(Ring, 0)], self.x, self.spaces, self.increasing), self

    def __mod__(self: Polynomial[Ring], other: Polynomial[Ring]) -> Polynomial[Ring]:
        """Return the remainder when dividing self by other.

        Examples:

            >>> numerator = Polynomial((1, 0, 0, 0, -1))
            >>> divisor = Polynomial((1, -1))
            >>> print(numerator); print(divisor)
            1 - x^4
            1 - x
            >>> print(numerator % divisor)
            0

            >>> numerator % Polynomial((1, 5))
            Traceback (most recent call last):
            ValueError: divisor must have leading coefficient 1 or -1, not 5

            >>> print(Polynomial((0,)) % divisor)
            0
        """
        #if not isinstance(other, Polynomial) and other.x_unwrapped == self.x_unwrapped:
        if not isinstance(other, Polynomial): #and other.x_unwrapped == self.x_unwrapped):
            return NotImplemented

        otherdeg = other._degree

        if otherdeg is None:
            raise ValueError("cannot divide by zero")

        if not isinstance(other, FPolynomial):
            if not (other[-1] == 1 or other[-1] == -1):
                raise ValueError(
                    f"divisor must have leading coefficient 1 or -1, not {other[-1]}"
                )

        #print("self", self)
        #print("other", other)
        #print("here1",self[-1] - self[-1] * other[-1], end=' ')
        #print("here2",self[-1] + self[-1] * other[-1], end=' ')
        #print("here3",(self[-1]-self[-1]*other[-1])*(self[-1]+self[-1]*other[-1]),end=' ')
        #print("type", type((self[-1]-self[-1]*other[-1])*(self[-1]+self[-1]*other[-1])))

        # debug some stuff
        factor1 = self[-1] - self[-1] * other[-1]
        factor2 = self[-1] + self[-1] * other[-1]
        has_to_be_zero = factor1 * factor2
        if isinstance(has_to_be_zero, Polynomial):
            assert has_to_be_zero._degree is None, (
                f"potential infinite loop in long division in divmod neither "
                f"{factor1} nor {factor2} are the zero poly, in Polynmial.mod()"
            )
        else:
            assert has_to_be_zero == cast(Ring, 0) * self[-1], (
                f"potential infinite loop in long division in divmod neither "
                f"{factor1} nor {factor2} are the zero poly, in Polynmial.mod()"
            )

        if self._degree is None:
            return self.__class__([cast(Ring, 0)], self.x, self.spaces, self.increasing)

        num = copy.copy(self)
        numdeg = num._degree

        if cast(int, numdeg) >= cast(int, otherdeg):
            while numdeg is not None and cast(int, numdeg) >= cast(int, otherdeg):
                monomial = Polynomial(
                    (cast(int, numdeg) - cast(int, otherdeg))
                    * [cast(Ring, 0)]
                    + [num[-1] * other[-1]], self.x, self.spaces, self.increasing
                )
                num = cast(
                    Polynomial, num - monomial * other
                )  # shouldn't have to cast here?
                numdeg = num._degree
            return num
        else:
            return self

    def __floordiv__(
        self: Polynomial[Ring], other: Polynomial[Ring]
    ) -> Polynomial[Ring]:
        """Return only the quotient when dividing self by other.

        Examples:

            >>> x = Polynomial([0, 1])
            >>> print((x**2 - 1) // (x - 1))
            1 + x

            >>> print((x - 1) // (x + 1 - x))
            -1 + x

            >>> from fractions import Fraction as F
            >>> x = FPolynomial([0, F(1)])
            >>> isinstance((x**2 - 1) // (x - 1), FPolynomial)
            True
        """
        return self.divmod(other)[0]

    def __hash__(self: Polynomial[Ring]) -> int:
        return hash(self._coeffs)

    def formalinv( self: Polynomial[Ring], maxdegree: int) -> Polynomial[Ring]:
        """Return formal (as series) inverse of self modulo  maxdegree+1.

        >>> x = Polynomial([0,1])
        >>> print((1 - x).formalinv(5))
        1 + x + x^2 + x^3 + x^4 + x^5
        >>> print((1 + x).formalinv(2))
        1 - x + x^2
        >>> print((1 - x**2).formalinv(5))
        1 + x^2 + x^4
        >>> print((1 - x**3).formalinv(6))
        1 + x^3 + x^6
        >>> print((-1 + x**2).formalinv(1))
        -1
        >>> print((x**12 - 1) * (x**3 - 1).formalinv(4))
        1 + x^3 - x^12 - x^15

        >>> p1 = (x**12 - 1)
        >>> p2 = (x**3 - 1)
        >>> deg_quot = p1._degree - p2._degree
        >>> (p1*p2.formalinv(deg_quot)).truncate(deg_quot) == p1 // p2
        True
        """
        if not (self[0] == 1 or self[0] == -1):
            raise ValueError(f"constant term must be 1 or -1, not {self[0]}")
        if not maxdegree > 0:
            raise ValueError(f"maxdegree should be positive, not {maxdegree}")

        firstnonzero = 1
        while self[firstnonzero] == 0:
            firstnonzero += 1
        realmax = maxdegree // firstnonzero

        accum = 1
        for _ in range(realmax):
            accum = (accum * (1 - self[0]*self) + 1).truncate(maxdegree)

        return self[0]*accum

    def derivative( self: Polynomial[Ring]) -> Polynomial[Ring]:
        """Return formal differentive of self.

        Example:

        >>> f = Polynomial([1,2,3,4])
        >>> print(f)
        1 + 2x + 3x^2 + 4x^3
        >>> print(f.derivative())
        2 + 6x + 12x^2
        """
        new_coeffs = []
        for i in range(self.degree()):
            new_coeffs.append((i+1)*self._coeffs[i+1])
        return self.__class__(new_coeffs, self.x, self.spaces, self.increasing)

    def apply(self, function):
        """Return copy of self with coefficient mapped according to function."""

        return self.__class__(tuple(map(function, self._coeffs)), self.x, self.spaces, self.increasing)

    def truncate(self, degree):
        """Return copy of self truncated beyond degree.

        Examples:

            >>> print(Polynomial([1]*5))
            1 + x + x^2 + x^3 + x^4
            >>> print(Polynomial([1]*5).truncate(2))
            1 + x + x^2
        """
        if degree < 0:
            raise ValueError(
                 f"truncation degree must be nonnegative, not {degree}"
            )

        return self.__class__(self._coeffs[:degree + 1], self.x, self.spaces, self.increasing)


Field = TypeVar("Field", bound=DivisionRing_)


class FPolynomial(Polynomial[Field]):

    """Extends the Polynomial class to polynomials over a field.

    The difference between the Polynomial class and this class is that,
    since we are working over a field, the denominator of the divisor
    in the divmod and __mod__ methods does not need to be monic.
    """
    __slots__ = ()

    def __init__(
        self: FPolynomial[Field], coeffs: Sequence[Field], x: str = "x", spaces: bool = True, increasing: bool = True
    ) -> None:
        """Create an polynomial over a field.

        The constructor here works the same as that of Polynomial. The difference
        between this class and the Polynomial class is that the divmod, __mod__,
        and __floordiv__ methods do not require other to be monic.

        Examples:

            >>> from fractions import Fraction
            >>> x = FPolynomial([0, Fraction(1)])
            >>> p1 = Fraction(-1,4) + Fraction(2, 3)*x + Fraction(-3, 5)*x**2
            >>> p2 = Fraction(7)*x + Fraction(1, 9)*x**2
            >>> print(p1 * p2)
            -7/4x + 167/36x^2 - 557/135x^3 - 1/15x^4
            >>> isinstance(p1 * p2, FPolynomial)
            True
            >>> isinstance(p1 // p2, FPolynomial)
            True
        """
        super().__init__(coeffs, x, spaces, increasing)

    def divmod(
        self: FPolynomial[Field], other: FPolynomial[Field]
    ) -> Tuple[FPolynomial[Field], ...]:

        """Return the quotient and remainder when dividing self by other.

        Examples:

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
        if not isinstance(other, FPolynomial):
            return NotImplemented

        if other._degree is None:
            raise ValueError("cannot divide by zero")

        if self._degree is None:
            return (
                self.__class__([cast(Field, 0)], self.x, self.spaces),
                self.__class__([cast(Field, 0)], self.x, self.spaces),
            )

        lead=other[-1]
        leadinv=lead**-1
        self_ = self.__class__(tuple(coef * leadinv for coef in self._coeffs), self.x, self.spaces, self.increasing)
        other_ = other.__class__(tuple(coef * leadinv for coef in other._coeffs), other.x, other.spaces, other.increasing)
        tup = super(FPolynomial, self_).divmod(other_)
        #tup = super(FPolynomial, self_ * other[-1] ** (-1)).divmod(other * other[-1] ** (-1))

        # should need to cast on next line?
        return (
            cast(FPolynomial[Field], tup[0]),
            #cast(FPolynomial[Field], tup[1] * other[-1]),
            cast(FPolynomial[Field], self.__class__(tuple(coef*lead for coef in tup[1]._coeffs),tup[1].x,tup[1].spaces,tup[1].increasing)),
        )

    def __mod__(
        self: FPolynomial[Field], other: FPolynomial[Field]
    ) -> FPolynomial[Field]:

        """Return the remainder when dividing self by other.

        Examples:

            >>> from fractions import Fraction
            >>> numerator = FPolynomial((1, 0, 0, 0, Fraction(-1)))
            >>> divisor = FPolynomial((1, Fraction(-1)))
            >>> print(numerator); print(divisor)
            1 - x^4
            1 - x
            >>> print(numerator % divisor)
            0

            >>> from fractions import Fraction
            >>> numerator = FPolynomial((1,0,0,0,Fraction(-1)))
            >>> divisor = FPolynomial((1, Fraction(-2)))
            >>> print(numerator); print(divisor)
            1 - x^4
            1 - 2x
            >>> print(numerator % divisor)
            15/16
        """
        if not isinstance(other, FPolynomial):
            return NotImplemented

        if other._degree is None:
            raise ValueError("cannot divide by zero")

        if self._degree is None:
            return cast(FPolynomial[Field], self.__class__([cast(Field, 0)], self.x, self.spaces))

        lead=other[-1]
        leadinv=lead**-1
        self_ = self.__class__(tuple(coef * leadinv for coef in self._coeffs), self.x, self.spaces, self.increasing)
        other_ = other.__class__(tuple(coef * leadinv for coef in other._coeffs), other.x, other.spaces, other.increasing)
        return cast(
            FPolynomial[Field],
            self.__class__(tuple(coef * lead for coef in super(FPolynomial, self_).__mod__(other_)), self.x, self.spaces, self.increasing)
        )

    def __repr__(self: FPolynomial[Field]) -> str:
        """Return repr string.

        Examples:

            >>> FPolynomial([])
            Traceback (most recent call last):
            ValueError: coeffs cannot be empty

            >>> FPolynomial([1, 2, 3, 0])
            FPolynomial((1, 2, 3))
        """
        return "FPolynomial(%s)" % repr(self._coeffs)

    # below here (and most if not all of the casting above is just for mypy purposes

    #def __add__(
    #    self: FPolynomial[Field], other: Union[Field, Polynomial[Field]]
    #) -> FPolynomial[Field]:
    #    return cast(FPolynomial[Field], super(FPolynomial, self).__add__(other))

    #def __radd__(self: FPolynomial[Field], other: Field) -> FPolynomial[Field]:
    #    if isinstance(other, DivisionRing_):
    #        return self.__class__((other,), self.x, self.spaces).__add__(self)
    #    else:
    #        return NotImplemented

    #def __neg__(self: FPolynomial[Field]) -> FPolynomial[Field]:
    #    return cast(FPolynomial[Field], super(FPolynomial, self).__neg__())

    #def __sub__(
    #    self: FPolynomial[Field], other: Polynomial[Field]
    #) -> FPolynomial[Field]:
    #    return cast(FPolynomial[Field], super(FPolynomial, self).__sub__(other))

    #def __rsub__(self: FPolynomial[Field], other: Field) -> FPolynomial[Field]:
    #    if isinstance(other, DivisionRing_):
    #        return self.__class__((other,), self.x, self.spaces).__sub__(self)
    #    else:
    #        return NotImplemented

    #def __mul__(
    #    self: FPolynomial[Field], other: Union[Field, Polynomial[Field]]
    #) -> FPolynomial[Field]:
    #    return cast(FPolynomial[Field], super(FPolynomial, self).__mul__(other))

    ## just calling super's rmul here causes some test to loop infinitely so
    ## redid all of the reverse ops above this way
    #def __rmul__(self: FPolynomial[Field], other: Field) -> FPolynomial[Field]:
    #    if isinstance(other, DivisionRing_):
    #        return self.__class__((other,), self.x, self.spaces).__mul__(self)
    #    else:
    #        return NotImplemented

    #def __pow__(self: FPolynomial[Field], n: int) -> FPolynomial[Field]:
    #    return cast(FPolynomial[Field], super(FPolynomial, self).__pow__(n))

    #def __floordiv__(
    #    self: FPolynomial[Field], other: Polynomial[Field]
    #) -> FPolynomial[Field]:
    #    return cast(FPolynomial[Field], super(FPolynomial, self).__floordiv__(other))

if __name__ == "__main__":

    import doctest

    doctest.testmod()
