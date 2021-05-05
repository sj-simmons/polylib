
import sys, unittest, copy

sys.path.insert(0, '..')
from polynomial import Polynomial, FPolynomial
from fractions import Fraction

class PolynomialTest(unittest.TestCase):

    def testZero(self):
        p = Polynomial([2,1,-2.3])
        self.assertEqual(str(p-p),'0')
        self.assertEqual((p-p).degree(),None)

    def testZero1(self):
        p = Polynomial([])
        self.assertEqual(str(p),'0')
        self.assertEqual(p.degree(),None)

    def testZero2(self):
        p = Polynomial([0])
        self.assertEqual(str(p),'0')
        self.assertEqual(p.degree(),None)

    def testZero3(self):
        p = Polynomial([.0])
        self.assertEqual(str(p),'0')
        self.assertEqual(p.degree(),None)

    def testZero4(self):
        p = Polynomial([3])
        self.assertEqual(str(p),'3')
        self.assertEqual(p.degree(),0)

    def testZero5(self):
        p = Polynomial([0])
        self.assertEqual(str(p-p),'0')
        self.assertEqual(p.degree(),None)

    def testZero6(self):
        p = Polynomial([])
        p = Polynomial([0])
        self.assertEqual(str(p-p),'0')
        self.assertEqual(p.degree(),None)

    def testContructor(self):
        p = Polynomial([2,1,-2.3])
        self.assertEqual(str(p),'2 + x - 2.3x^2')
        self.assertEqual(p.degree(),2)

    def testContructor2(self):
        p = Polynomial([0,-3,2.3,1/4])
        self.assertEqual(str(p),'-3x + 2.3x^2 + 0.25x^3')
        self.assertEqual(p.degree(),3)

    def testContructorDefault(self):
        p = Polynomial([])
        self.assertEqual(str(p),'0')
        self.assertEqual(p.degree(),None)

    def testContructorDefault2(self):
        p = Polynomial(())
        self.assertEqual(str(p),'0')
        self.assertEqual(p.degree(),None)

    def testAddition(self):
        p1 = Polynomial([2,1,3])
        p2 = Polynomial([0,-1,-.5,0,4])
        self.assertEqual(str(p1+p2),'2 + 2.5x^2 + 4x^4')
        self.assertEqual((p1+p2).degree(),4)

    def testAdditionZero(self):
        p1 = Polynomial([2,1,3])
        p2 = Polynomial([0])
        self.assertEqual(str(p1+p2),'2 + x + 3x^2')
        self.assertEqual((p1+p2).degree(),2)

    def testDegree(self):
        p = Polynomial([2,3,0])
        self.assertEqual(p.degree(),1)

    def testDegree2(self):
        p1 = Polynomial([2,3,2,-1])
        p2 = Polynomial([2,3,-2,1])
        self.assertEqual(str(p1+p2),'4 + 6x')
        self.assertEqual((p1+p2).degree(),1)

    def testOnes(self):
        p = Polynomial([0,-1,1.0])
        self.assertEqual(str(p),'-x + x^2')

    def testSubtraction(self):
        p1 = Polynomial([2,1,3,0,4])
        p2 = Polynomial([0,-1,-0.5,0.23,4])
        self.assertEqual(str(p1-p2),'2 + 2x + 3.5x^2 - 0.23x^3')
        self.assertEqual((p1-p2).degree(),3)
        self.assertEqual(p2,p2)   # make sure negation didn't change p2

    def testSubtraction2(self):
        p1 = Polynomial([2, 1, 3, 0, 4.0])  #2 + x + 3x^2 + 4x^4
        p2 = Polynomial([0,-1,-0.5,0.23,4.0])  #-x - .5 x^2 + .23 x^3 + 4x^4
        self.assertEqual(str(p1-p2),'2 + 2x + 3.5x^2 - 0.23x^3')
        self.assertEqual((p1-p2).degree(),3)
        self.assertEqual(p2,p2)   # make sure negation didn't change p2

    def testMultiplication(self):
        p1 = Polynomial([2, 1, 3]) # 2 + x + 3x^2
        p2 = Polynomial([0, 1, 0])   # x
        self.assertEqual(str(p1*p2),'2x + x^2 + 3x^3')
        self.assertEqual((p1*p2).degree(),3)

    def testMultiplication2(self):
        p1 = Polynomial([2, 1, 3])  # 2 + x + 3x^2
        p2 = Polynomial([0,-1, 4])  # -x + 4x^2
        self.assertEqual(str(p1*p2),'-2x + 7x^2 + x^3 + 12x^4')
        self.assertEqual((p1*p2).degree(),4)

    def testMultiplication3(self):
        p = Polynomial([1, -1])
        self.assertEqual(str(p*p),'1 - 2x + x^2')

    def testScalarMultiplication4(self):
        p = Polynomial([1, -1, .3])
        self.assertEqual(str(2*p),'2 - 2x + 0.6x^2')

    def testScalarMultiplication5(self):
        scalar = Polynomial([.2])    #can you have real scalar mult by non-ints?
        p = Polynomial([1, -1, .3])
        self.assertEqual(str(scalar*p),'0.2 - 0.2x + 0.06x^2')

    def testScalarMultiplication6(self):
        scalar = Polynomial([.2])    #can you have real scalar mult by non-ints?
        p = Polynomial([1, -1, .3])
        self.assertEqual(str(p*scalar),'0.2 - 0.2x + 0.06x^2')

    def testScalarAdd(self):
        scalar = Polynomial([.2])
        p = Polynomial([1, -1, .3])
        self.assertEqual(str(scalar + p),'1.2 - x + 0.3x^2')

    def testScalarAdd2(self):
        scalar = Polynomial([.2])
        p = Polynomial([1, -1, .3])
        self.assertEqual(str(p + scalar),'1.2 - x + 0.3x^2')

    def testFractionConstructor(self):
        p = Polynomial([Fraction(1), Fraction(-1,3)])
        self.assertEqual(str(p),'1 - 1/3x')

    def testFractionToZero(self):
        p = Polynomial([Fraction(1), Fraction(-1,3)])
        self.assertEqual(str(p**0),'1')

    def testFractionMulitplication(self):
        p1 = Polynomial([Fraction(1), Fraction(-1,3)])
        p2 = Polynomial([2])
        self.assertEqual(str(p1*p2),'2 - 2/3x')

    def testFractionMulitplication2(self):
        p1 = Polynomial([Fraction(1), Fraction(-1,3)])
        p2 = Polynomial([2])
        self.assertEqual(str(p2*p1),'2 - 2/3x')

    def testScalarPower(self):
        p = Polynomial([Fraction(2)])
        self.assertEqual(str(p**3),'8')

    def testEvaluation(self):
        p = Polynomial([1,Fraction(1,2)])
        self.assertEqual(str(p.of(3)),'5/2')

    def testEvaluation2(self):
        p = Polynomial([1,Fraction(1,2)])
        self.assertEqual(str(p.of(Fraction(3,2))),'7/4')

    def testCopy(self):
        p = Polynomial([1,2,3])
        pnew = copy.copy(p)
        p._coeffs = Polynomial([10,11,12]) # but don't ever do this
        self.assertEqual(str(pnew),'1 + 2x + 3x^2')

class FPolynomialTest(unittest.TestCase):

    def testType(self):
        p = FPolynomial([2,1,3])
        self.assertEqual(str(type(p)),"<class 'polynomial.FPolynomial'>")

    def testType2(self):
        p = FPolynomial([2,1,3])
        pnew = copy.copy(p)
        self.assertEqual(str(type(pnew)),"<class 'polynomial.FPolynomial'>")

    def testType3(self):
        p1 = FPolynomial([2,1,3])
        p2 = FPolynomial([1,3])
        sum = p1 + p2
        self.assertEqual(str(type(sum)),"<class 'polynomial.FPolynomial'>")

    def testAddition(self):
        p1 = FPolynomial([2,1,3])
        p2 = FPolynomial([0,-1,-.5,0,4])
        self.assertEqual(str(p1+p2),'2 + 2.5x^2 + 4x^4')
        self.assertEqual((p1+p2).degree(),4)

    def testDivmod1(self):
        p = FPolynomial([1,2,3])
        q, r = p.divmod(2)
        self.assertEqual(str(q),'0.5 + x + 1.5x^2')
        self.assertEqual(str(r),'0')

    def testDivmod2(self):
        p = FPolynomial([Fraction(1),Fraction(2),Fraction(3)])
        q, r = p.divmod(2)
        self.assertEqual(str(q),'1/2 + x + 3/2x^2')
        self.assertEqual(str(r),'0')

    def testDivmod2(self):
        p = FPolynomial([1,0,0,0,-1])
        q, r = p.divmod(FPolynomial([1,-1]))
        self.assertEqual(str(q),'1.0 + x + x^2 + x^3')
        self.assertEqual(str(r),'0')

    #Fails due to rounding error
    #def testDivmod2(self):
    #    p1 = FPolynomial([1,2.1,7,-2,-1])
    #    p2 = FPolynomial([1,-1,3.2])
    #    q, r = p1.divmod(p2)
    #    self.assertEqual(p1, q * p2 + r)

    def testDivmod2(self):
        p1 = FPolynomial([Fraction(1),Fraction(1,2),Fraction(7),-2,-1])
        p2 = FPolynomial([1,-1,Fraction(-3,4)])
        q, r = p1.divmod(p2)
        self.assertEqual(p1, q * p2 + r)

def main(argv):
    unittest.main()

if __name__ == '__main__':
    main(sys.argv)
