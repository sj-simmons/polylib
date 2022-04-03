from polylib import Polynomial
from functools import reduce
from operator import mul

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


def cyclotomic(n: int, moebius=False, gmp=False) -> Polynomial:
    """Return the nth cyclotomic polynomial.

    For a positive integer n, this returns the nth cyclotomic polynomial
    with coefficients of type int (or of thype mpz if gmp is True).

    If numlib is installed, this computes via the Moebius inversion form-
    ula, which can be much faster.

    Note: gmp mpzs are faster than Python's longs "once the precision ex-
    ceeds 20 to 50 bits."  Unless you know that the coefficients are that
    big, then leave gmp set to False.

    Args:

        n (int): the cyclotomic polynomial to return.

        moebius (bool): if True, use the faster method to compute the cyc-
            lotomic polynomial. Requires numlib being installed.

        gmp: If True, use mpz integers instead of Python longs. Requires
            gmppy2 being installed.

    Returns:

        (Polynomial). The cyclotomic polynomial with coefficients of type
            Python int.



    Examples:

        >>> print(cyclotomic(60))
        1 + x^2 - x^6 - x^8 - x^10 + x^14 + x^16

        The polynomial cyclotomic(n) has integer coeffients for all n. There
        are multiple ways to, say, reduce those coefficients modulo a prime or,
        more generally, to map to the coefficients to come other ring of field.

        If you want coefficients but of type Fraction:

        >>> from polylib import FPolynomial
        >>> from fractions import Fraction
        >>> x = FPolynomial([0, Fraction(1)]) # an indeterminant for Q[x]
        >>> print(cyclotomic(2)(x))
        1 + x
        >>> print(repr(cyclotomic(2)(x)))
        FPolynomial((Fraction(1, 1), Fraction(1, 1)))

        If you want to reduce the coefficients modulo n:

        >>> from polylib import Polynomial
        >>> from numlib import Zmod
        >>> GF = Zmod(17) # Z/17Z
        >>> cyclotomic(2) * GF(1)
        Polynomial((1 + <17>, 1 + <17>))

        In the last example, the coefficients are now in a field; so you
        may wish to recast to FPolynomial

        >>> FPolynomial(cyclotomic(2) * GF(1))
        FPolynomial((1 + <17>, 1 + <17>))

        But you may be better off, speed-wise, to first compute in the
        integers and then specialize the coefficients.
        >>> poly = cyclotomic(1001)
        >>> poly  # doctest: +ELLIPSIS
        Polynomial((1, 1, 1, ...
        >>> #poly.apply(lambda coeff: Zmod(17)(coeff))
        >>> poly * Zmod(17)(1)  # doctest: +ELLIPSIS
        Polynomial((1 + <17>, 1 + <17>, 1 + <17>, ...
        >>> FPolynomial(poly * Zmod(17)(1))  # doctest: +ELLIPSIS
        FPolynomial((1 + <17>, 1 + <17>, 1 + <17>, ...
        >>> #cyclopoly = cyclotomic(2**12-1)
        >>> #cyclopoly.apply(lambda coeff: Zmod(17)(coeff))
        >>> #cyclopoly.degree()
        1728
    """
    if not isinstance(n, int):
        raise TypeError(f"n must be of type int, not {type(n)}")
    if n < 1:
        raise ValueError(f"n must be positive, not {n}")

    if moebius:
        # try:
        # # below computes literally the nth cyclotomic poly via its
        # # moebius inversion formula.  In practice, this is slower
        # # than the recursive formula in general, and particularly
        # # on products of say 4 or more distinct primes.
        # from numlib import mu
        # prod = Polynomial([-1] + [0]*(n-1) + [1])
        # for d in range(1, n):
        #     if n % d == 0:
        #         if mu(n//d) == 1:
        #             prod *= Polynomial([-1] + [0]*(d-1) + [1])
        #         elif mu(n//d) == -1:
        #             prod //= Polynomial([-1] + [0]*(d-1) + [1])
        # return prod

        import numlib  # type: ignore

        def cyclo_distinct_primes(n, one=1):
            """Return cyclo. poly. for n = prod. of distinct primes."""
            # below computes literally the nth cyclotomic poly via its
            # moebius inversion formula.  In practice, this is slower
            # than the recursive formula in general, and particularly
            # on products of say 4 or more distinct primes.
            from numlib import mu

            deg = numlib.phi(n)

            # prod = 1
            # for d in range(1, n):
            #    if n % d == 0:
            #        if numlib.mu(n//d) == 1:
            #            prod *= Polynomial([-1] + [0]*(d-1) + [1])
            #        elif numlib.mu(n//d) == -1:
            #            prod*(Polynomial([-1] + [0]*(d-1) + [1])).formalinv(deg)
            # return prod

            nums = []
            dens = []
            for d in range(1, n + 1):
                if n % d == 0:
                    if numlib.mu(n // d) == 1:
                        nums.append(Polynomial([one] + [0 * one] * (d - 1) + [one]))
                    else:
                        dens.append(Polynomial([one] + [0 * one] * (d - 1) + [one]))
            # print(', '.join(map(lambda x: str(x), nums)))
            # print(', '.join(map(lambda x: str(x), dens)))

            # slow
            # num = reduce(mul, nums, 1)
            # den = reduce(mul, dens, 1)
            # return (num.truncate(deg) * (den.formalinv(deg))).truncate(deg)

            # still slow
            num = reduce(mul, nums, 1)
            return reduce(
                lambda x, y: (x * (y.formalinv(x._degree - y._degree))).truncate(
                    x._degree - y._degree
                ),
                dens,
                num,
            )

            # fftmult not faster, here, for some reason (poly mult implementation
            # might already be nlog(n)??)
            # on 5*7*11*13,  nearly 51.1 sec (per loop)
            # on 3*5*7*11*13,
            # return reduce(
            #    lambda x,y: (x*y).truncate(deg), [
            #        num.fftmult(den.formalinv(deg)).truncate(deg)
            #        for num, den in zip(nums, dens)
            #    ], 1
            # ).apply(round).apply(int)

            #                   Core i5 (int)
            # on 5*7*11*13,     5.88 secs (per loop)
            # on 3*5*7*11*13,   54   secs (per loop)
            # on 2*3*5*7*11*13     secs (per loop)
            # on 255255              mins
            return reduce(
                lambda x, y: (x * y).truncate(deg),
                [
                    (num * den.formalinv(deg)).truncate(deg)
                    for num, den in zip(nums, dens)
                ],
                one,
            )

            # on 5*7*11*13,      8.45  sec (per loop)
            # on 3*5*7*11*13,    1m16s
            # on 2*3*5*7*11*13   3m31s   WHY IS THIS SLOWER THAN ABOVE
            # on 255255
            # consider truncating during multiplication in both places below; so
            # add truncate paramter to __mul__ etc.
            # return reduce(
            #    lambda x,y: (x*y).truncate(deg), [
            #        (num * den.formalinv(num._degree-den._degree)).truncate(num._degree-den._degree)
            #        if num._degree % den._degree == 0 else
            #        (num * den.formalinv(deg)).truncate(deg)
            #        for num, den in zip(nums, dens)
            #    ], 1
            # )

            #                   Core i5 (int)
            # on 5*7*11*13,     6.23 secs (per loop)
            # on 3*5*7*11*13,   53.6 secs (per loop)
            # on 2*3*5*7*11*13  141  secs (per loop)
            # on 255255         460  mins
            # return reduce(
            #    lambda x,y: (x*y).truncate(deg), [
            #        (num * den.formalinv(num._degree-den._degree)).truncate(num._degree-den._degree)
            #        if num._degree % den._degree == 0 else
            #        (num * den.formalinv(deg - (deg % den._degree))).truncate(deg)
            #        if den._degree <= deg else
            #        num.truncate(deg)
            #        for num, den in zip(nums, dens)
            #    ], one
            # )

        one = 1

        if gmp:
            try:
                import gmpy2  # type: ignore

                one = gmpy2.mpz(1)
            except:
                pass

        return cyclo_distinct_primes(n, one).apply(int)

    # except:
    #    pass

    # on 5*7*11*13, 13.75 mins
    # there is no point in putting option to use gmp here
    def cyclopoly(n: int) -> Polynomial:
        return (Polynomial([-1] + [0] * (n - 1) + [1])) // reduce(
            mul, [cyclopoly(d) for d in range(1, n) if n % d == 0], Polynomial([1])
        )

    return cyclopoly(n)


if __name__ == "__main__":

    import doctest

    doctest.testmod()
