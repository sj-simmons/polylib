from __future__ import annotations
import math
import numlib as nl
#from numlib import mu, phi, factor

def intpolypow(p: IntPoly, n: int) -> IntPoly:

    if n == 0:
        return IntPoly((1,))

    fact = intpolypow(p, n // 2)
    if n % 2 == 0:
        return fact * fact

    return fact * fact * p

class IntPoly:
    def __init__(self, coeffs: tuple[int, ...]):

        while coeffs[-1] == 0:
            coeffs = coeffs[:-1]
        self.cos = coeffs

    def __add__(self, __x: IntPoly) -> IntPoly:

        selflen = len(self.cos)
        xlen = len(__x.cos)

        if selflen >= xlen:
            return IntPoly(
                tuple(self.cos[i] + __x.cos[i] for i in range(xlen))
                + self.cos[xlen:]
            )
        else:
            return IntPoly(
                tuple(self.cos[i] + __x.cos[i] for i in range(selflen))
                + __x.cos[selflen:]
            )

    def __sub__(self, __x: IntPoly) -> IntPoly:

        selflen = len(self.cos)
        xlen = len(__x.cos)

        if selflen >= xlen:
            return IntPoly(
                tuple(self.cos[i] - __x.cos[i] for i in range(xlen))
                + self.cos[xlen:]
            )
        else:
            return IntPoly(
                tuple(self.cos[i] - __x.cos[i] for i in range(selflen))
                + tuple(-y for y in __x.cos[selflen:])
            )

    def __neg__(self) -> IntPoly:

        return IntPoly(tuple(-y for y in self.cos))

    def __mul__(self, __x: IntPoly) -> IntPoly:

        selflen = len(self.cos)
        xlen = len(__x.cos)

        product = []

        if selflen >= xlen:
            for i in range(selflen):
                sum_ = 0
                if i < xlen:
                    for j in range(i + 1):
                        sum_ += __x.cos[j] * self.cos[i - j]
                    product.append(sum_)
                else:
                    for j in range(xlen):
                        sum_ += __x.cos[j] * self.cos[i - j]
                    product.append(sum_)
            for i in range(xlen - 1):
                sum_ = 0
                for j in range(i + 1, xlen):
                    sum_ += __x.cos[j] * self.cos[selflen + i - j]
                product.append(sum_)
        else:
            for i in range(xlen):
                sum_ = 0
                if i < selflen:
                    for j in range(i + 1):
                        sum_ += self.cos[j] * __x.cos[i - j]
                    product.append(sum_)
                else:
                    for j in range(selflen):
                        sum_ += self.cos[j] * __x.cos[i - j]
                    product.append(sum_)
            for i in range(selflen - 1):
                sum_ = 0
                for j in range(i + 1, selflen):
                    sum_ += self.cos[j] * __x.cos[xlen + i - j]
                product.append(sum_)

        return IntPoly(tuple(product))

    def __xor__(self, __n: int) -> IntPoly:

        return intpolypow(self, __n)

    def __str__(self) -> str:
        var = "x"
        s = ""
        for i in range(len(self.cos)):
            coef = self.cos[i]
            if coef > 0:  # add coefficient
                if i != 0 and coef == 1 and s != "":
                    s += " + "
                elif i != 0 and s != "":
                    s += " + " + str(coef)
                elif i == 0 or coef != 1:
                    s += str(coef)
            elif coef < 0:
                if coef == -1 and s == "":
                    if i == 0:
                        s += "-1"
                    else:
                        s += "-"
                elif coef == -1:
                    s += " - "
                elif s == "":
                    s += "-" + str(-coef)
                else:
                    s += " - " + str(-coef)
            if i > 1 and coef != 0:  # add x^n
                s += var + "^" + str(i)
            elif i == 1 and coef != 0:  # add x
                s += var
        return s

    def __repr__(self) -> str:
        return str(self)

    def formalinv(self, maxdegree: int) -> IntPoly:
        """Return truncated inverse 1/self of self with constant term 1 or -1."""

        const = self.cos[0]
        cos = self.cos + (0,) * (maxdegree - len(self.cos) + 1)
        inv = [const]

        for n in range(1, maxdegree + 1):
            accum = 0
            for k in range(n):
                accum += -inv[k] * const * cos[n-k]
            inv.append(accum)

        return IntPoly(tuple(inv))

    def __call__(self, x: int) -> int:
        b = self.cos[-1]
        for co in self.cos[:-1:-1]:
            b = co + b * x
        return b

    #def divmod(self, __x: IntPoly) -> tuple(IntPoly, IntPoly):
    #    """Return quotient and remainder when self by monic __x."""

    #    if self == IntPoly((0,)):
    #        return (IntPoly((0,)), IntPoly((0,)))

    #    while num[-1] == 0:
    #        num  = num[:-1]

    #    num = self.cos
    #    quo = (0,) * (len(num.cos) - len(__x.cos) - 1) + num.cos[-1]

    #    num -=  quo * __x
    #    while len(rem.cos) >= len(num.cos):
    #        num = self + IntPoly((-self[0]),) * __x

    #def __floordiv__(self, __x: IntPoly) -> IntPoly:
    #    """Return quotient after dividing self by monic __x using div. algo."""
    #    #pass



class Rational:
    """
    Rational number class.

    This class compiles but would performs no better than fractions.Fraction
    were it not for the fact that it uses non-reduced represations.  The re-
    sult is far less calls to math.gcd for programs with tens of thousands of
    rational number operations.

    One must remember that non-reduced representations are being in play and,
    when necessary, reduce using the reduce() method.  Though string and repr
    are automatically pre-reduced.
    """

    def __init__(self, a: int, b: int = 1, bitshift: int = 3000):

        if bitshift and a >> bitshift > 0 or b >> bitshift > 0:
            gcd_ = math.gcd(a, b)
            self.a = a // gcd_
            self.b = b // gcd_
        else:
            self.a = a
            self.b = b

    def __add__(self, __x: Rational) -> Rational:

        return Rational(self.a * __x.b + __x.a * self.b, self.b * __x.b)

    def __sub__(self, __x: Rational) -> Rational:

        return Rational(self.a * __x.b - __x.a * self.b, self.b * __x.b)

    def __neg__(self) -> Rational:

        return Rational(-self.a, self.b)

    def __mul__(self, __x: Rational) -> Rational:

        return Rational(self.a * __x.a, self.b * __x.b)

    def __truediv__(self, __x: Rational) -> Rational:

        return Rational(self.a * __x.b, self.b * __x.a)

    def __xor__(self, __n: int) -> Rational:

        if __n >= 0:
            return Rational(self.a ** __n, self.b ** __n)
        else:
            return Rational(self.b ** -__n, self.a ** -__n)

    def __lt__(self, __x: Rational) -> bool:

        return self.a * __x.b < self.b * __x.a

    def __le__(self, __x: Rational) -> bool:

        return self.a * __x.b <= self.b * __x.a

    def __gt__(self, __x: Rational) -> bool:

        return self.a * __x.b > self.b * __x.a

    def __ge__(self, __x: Rational) -> bool:

        return self.a * __x.b >= self.b * __x.a

    def __eq__(self, __x: object) -> bool:

        if not isinstance(__x, Rational):
            return NotImplemented
        return self.a * __x.b == self.b * __x.a

    def __str__(self) -> str:

        gcd_ = math.gcd(self.a, self.b)
        self.a //= gcd_
        self.b //= gcd_

        if self.b > 1:
            return f"{self.a}/{self.b}"

        if self.b < -1:
            return f"{-self.a}/{-self.b}"

        return f"{self.b * self.a}"

    def __repr__(self) -> str:
        return str(self)

    def reduced(self):

        gcd_ = math.gcd(self.a, self.b)
        self.a //= gcd_
        self.b //= gcd_

        return self


def rationalpolypow(p: RationalPoly, n: int) -> RationalPoly:

    if n == 0:
        return RationalPoly((Rational(1),))

    fact = rationalpolypow(p, n // 2)
    if n % 2 == 0:
        return fact * fact

    return fact * fact * p


class RationalPoly:
    def __init__(self, coeffs: tuple[Rational, ...]):

        while coeffs[-1] == Rational(0):
            coeffs = coeffs[:-1]
        self.cos = coeffs

    def __add__(self, __x: RationalPoly) -> RationalPoly:

        selflen = len(self.cos)
        xlen = len(__x.cos)

        if selflen >= xlen:
            return RationalPoly(
                tuple(self.cos[i] + __x.cos[i] for i in range(xlen))
                + self.cos[xlen:]
            )
        else:
            return RationalPoly(
                tuple(self.cos[i] + __x.cos[i] for i in range(selflen))
                + __x.cos[selflen:]
            )

    def __sub__(self, __x: RationalPoly) -> RationalPoly:

        selflen = len(self.cos)
        xlen = len(__x.cos)

        if selflen >= xlen:
            return RationalPoly(
                tuple(self.cos[i] - __x.cos[i] for i in range(xlen))
                + self.cos[xlen:]
            )
        else:
            return RationalPoly(
                tuple(self.cos[i] - __x.cos[i] for i in range(selflen))
                + tuple(-y for y in __x.cos[selflen:])
            )

    def __neg__(self) -> RationalPoly:

        return RationalPoly(tuple(-y for y in self.cos))

    def __mul__(self, __x: RationalPoly) -> RationalPoly:

        selflen = len(self.cos)
        xlen = len(__x.cos)

        product = []

        if selflen >= xlen:
            for i in range(selflen):
                sum_ = Rational(0)
                if i < xlen:
                    for j in range(i + 1):
                        sum_ += __x.cos[j] * self.cos[i - j]
                    product.append(sum_)
                else:
                    for j in range(xlen):
                        sum_ += __x.cos[j] * self.cos[i - j]
                    product.append(sum_)
            for i in range(xlen - 1):
                sum_ = Rational(0)
                for j in range(i + 1, xlen):
                    sum_ += __x.cos[j] * self.cos[selflen + i - j]
                product.append(sum_)
        else:
            for i in range(xlen):
                sum_ = Rational(0)
                if i < selflen:
                    for j in range(i + 1):
                        sum_ += self.cos[j] * __x.cos[i - j]
                    product.append(sum_)
                else:
                    for j in range(selflen):
                        sum_ += self.cos[j] * __x.cos[i - j]
                    product.append(sum_)
            for i in range(selflen - 1):
                sum_ = Rational(0)
                for j in range(i + 1, selflen):
                    sum_ += self.cos[j] * __x.cos[xlen + i - j]
                product.append(sum_)

        return RationalPoly(tuple(product))

    def __xor__(self, __n: int) -> RationalPoly:

        return rationalpolypow(self, __n)

    def __str__(self) -> str:
        var = "x"
        s = ""
        for i in range(len(self.cos)):
            coef = self.cos[i]
            if coef > Rational(0):  # add coefficient
                if i != 0 and coef == Rational(1) and s != "":
                    s += " + "
                elif i != 0 and s != "":
                    s += " + " + str(coef)
                elif i == 0 or coef != Rational(1):
                    s += str(coef)
            elif coef < Rational(0):
                if coef == Rational(-1) and s == "":
                    if i == 0:
                        s += "-1"
                    else:
                        s += "-"
                elif coef == Rational(-1):
                    s += " - "
                elif s == "":
                    s += "-" + str(-coef)
                else:
                    s += " - " + str(-coef)
            if i > 1 and coef != Rational(0):  # add x^n
                s += var + "^" + str(i)
            elif i == 1 and coef != Rational(0):  # add x
                s += var
        return s

    def __repr__(self) -> str:
        return str(self)

    def formalinv(self, maxdegree: int) -> RationalPoly:
        """Return truncated inverse 1/self of self with nonzero constant term."""

        const = self.cos[0]
        cos = self.cos + (Rational(0),) * (maxdegree - len(self.cos) + 1)
        inv = [Rational(1)/const]

        #for n in range(1, maxdegree + 1):
        for n in range(1, maxdegree + 1):
            accum = Rational(0)
            for k in range(n):
                accum += -inv[k] / const * cos[n-k]
            inv.append(accum)

        return RationalPoly(tuple(inv))

    def reduced(self):

        return RationalPoly(tuple(c.reduced() for c in self.cos))

def bernoulli(n: int) -> RationalPoly:
    """
    Return the generating series x/(1-e^(-x)) modulo (x^(n+1)).

    One has: Sum_{i=0} B_n (-x)^n/n! = x/(1-e^(-x)).

    Interactively, e.g.,

    >>> print(bernoulli(10))
    1 + 1/2x + 1/12x^2 - 1/720x^4 + 1/30240x^6 - 1/1209600x^8 + 1/47900160x^10

    The corresponding Bernoulli numbers:

    >>> import math
    >>> print(", ".join(
    ...     f"B_{i} = {Rational((-1)**i*math.factorial(i))*x}"
    ...     for i, x in enumerate(bernoulli(10).cos))
    ... ) # doctest: +ELLIPSIS
    B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_3 = 0, B_4 = -1/30, ...
    """
    if n == 0:
        return RationalPoly((Rational(1),))

    # generate n terms of the series p(x) satisfying x/(1-e^(-x)) = 1/(1-p(x))
    p_ = [Rational(1)]
    for i in range(2, n + 2):
        p_.append(-Rational((-1) ** i, math.factorial(i)))

    return RationalPoly(tuple(p_)).formalinv(n).reduced()

#def fft(tuple[int]):
#    return 

class SparseIntPoly:

    def __init__(self, coeffs: tuple[int], sparseness: int = 0):
        self.cos = coeffs
        self.sparse = sparseness

    def __mul__(self, __x):
        pass

    def __truediv__(self, __x):
        pass

    def toIntPoly(self):
        pass


def oddsqfree(n: int) -> IntPoly:
    """Helper for cyclotomic.

    Returns the cyclotomic polynomial for an odd, squarefree integer
    greater than 1.
    """
    nums: list[IntPoly] = []
    dens: list[IntPoly] = []

    for d in range(1, n+1):
        if n % d == 0:
            if nl.mu(n//d) == 1:
                #nums.append(IntPoly((-1,) + (0,)*(d-1) + (1,)))
                nums.append(IntPoly((1,) + (0,)*(d-1) + (-1,)))
            elif nl.mu(n//d) == -1:
                #dens.append(IntPoly((-1,) + (0,)*(d-1) + (1,)))
                dens.append(IntPoly((1,) + (0,)*(d-1) + (-1,)))

    prod = IntPoly((1,))
    phi_ = nl.phi(n)
    half = phi_//2 + 1
    for num, den in zip(nums, dens):
        prod = IntPoly((prod * IntPoly((num * den.formalinv(phi_)).cos[:half])).cos[:half])

    cos = prod.cos + (0,) * (half - len(prod.cos))
    return IntPoly(cos + cos[-2::-1])
    #return IntPoly(prod.cos + prod.cos[-2::-1])

    #prod = IntPoly((1,))
    #for d in range(1, n+1):
    #    if n % d == 0:
    #        if mu(n//d) == 1:
    #            #prod *= IntPoly((-1,) + (0,)*(d-1) + (1,))
    #            if d > phi_:
    #                prod =  prod * IntPoly((-1,))
    #            else:
    #                prod =  prod * IntPoly((-1,)) + IntPoly(((0,) * d + prod.cos)[:phi_+1])
    #        elif mu(n//d) == -1:
    #            #prod //= IntPoly([-1] + [0]*(d-1) + [1])
    #            prod *= IntPoly((-1,) + (0,)*(d-1) + (1,)).formalinv(phi_)
    #            prod = IntPoly(prod.cos[:phi_+1])  #NOTE: can this be improved
    #return prod

    #if n == 1:
    #    return IntPoly((-1, 1))
    #prod = IntPoly((1,))
    #for d in range(1, n):
    #    prod *= IntPoly((-1,) + (0,)*(n-1) + (1,)) // cyclotomic(d)
    #return prod

def cyclotomic(n: int) -> IntPoly:

    if n == 1:
        return IntPoly((-1, 1))
    if n == 2:
        return IntPoly((1, 1))

    for p, e in nl.factor2(n):
        if e > 1:
            newcos = (1,)
            for co in cyclotomic(n//p).cos[1:]:
                newcos += (0,)*(p-1) + (co,)
            return IntPoly(newcos)
        if p == 2:
            newcos = ()
            for i, co in enumerate(cyclotomic(n//2).cos):
                newcos += (co,) if i % 2 == 0 else (-co,)
            return IntPoly(newcos)

    return oddsqfree(n)


if __name__ == '__main__':

    import doctest
    doctest.testmod()
