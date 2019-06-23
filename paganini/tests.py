import unittest

from paganini.expressions import *
from paganini.specification import *

class SingularTuner(unittest.TestCase):

    def test_singular_btrees(self):
        """ Singular tuning of binary trees
            B = 1 + Z * B^2."""

        spec = Specification()
        z, B = Variable(), Variable()
        spec.add(B, 1 + z * B ** 2)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.25)
        self.assertAlmostEqual(B.value, 2)

    def test_singular_motzkin_trees(self):
        """ Singular tuning of Motzkin trees
            M = Z * SEQ<=2(M). """

        spec = Specification()
        z, M = Variable(), Variable()
        spec.add(M, z * Seq(M, leq(2)))

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.333333333333334)
        self.assertAlmostEqual(M.value, 1.0)

    def test_singular_motzkin_trees2(self):
        """ Singular tuning of Motzkin trees
            M = Z + Z * M + Z * M^2. """

        spec = Specification()
        z, M = Variable(), Variable()
        spec.add(M, z + z * M + z * M ** 2)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.333333333333334)
        self.assertAlmostEqual(M.value, 1.0)

    def test_singular_trees(self):
        """ Singular tuning of plane trees
            T = Z * SEQ(T)."""

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, z * Seq(T))

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.25)
        self.assertAlmostEqual(T.value, 0.5)

    def test_singular_lambda_terms(self):
        """ Singular tuning of plain lambda terms
            L = Z * SEQ(Z) + Z * L + Z * L^2."""

        spec = Specification()
        z, L, D = Variable(), Variable(), Variable()
        spec.add(L, D + z * L + z * L ** 2)
        spec.add(D, z + z * D)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.295597742522085)
        self.assertAlmostEqual(L.value, 1.19148788395312)

    def test_singular_lambda_terms2(self):
        """ Singular tuning of plain lambda terms
            L = Z * SEQ(Z) + Z * L + Z * L^2."""

        spec = Specification()
        z, L = Variable(), Variable()
        spec.add(L, z * Seq(z) + z * L + z * L ** 2)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.295597742522085)
        self.assertAlmostEqual(L.value, 1.19148788395312)

    def test_singular_polya_trees(self):
        """ Singular tuning of Polya trees
            T = Z * MSet(T)."""

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, z * MSet(T))

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.338322112871298)
        self.assertAlmostEqual(T.value, 1)

    def test_singular_custom_trees(self):
        """ Singular tuning of some custom trees defined by
            T = Z + Z * SEQ_>=2(T)."""

        params = Params(Type.ALGEBRAIC)
        params.max_iters = 100 # required

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, z + z * Seq(T, geq(2)))

        spec.run_singular_tuner(z, params)

        self.assertAlmostEqual(z.value, 0.333333333333335)
        self.assertAlmostEqual(T.value, 0.499999999999993)

    def test_binary_words(self):
        """ Singular tuning of binary words.
            B = SEQ(Z + Z). """

        spec = Specification()
        z, B = Variable(), Variable()
        spec.add(B, Seq(z + z))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.5, 5)

    def test_compositions(self):
        """ Singular tuning of all compositions.
            C = SEQ(Z * SEQ(Z)). """

        params = Params(Type.RATIONAL)
        params.max_iters = 10000 # required

        spec = Specification()
        z, C = Variable(), Variable()
        spec.add(C, Seq(z * Seq(z)))

        spec.run_singular_tuner(z, params)
        self.assertAlmostEqual(z.value, 0.5, 4)

    def test_compositions_with_restricted_summands(self):
        """ Singular tuning of compositions with restricted summands in {1,2}.
            C = SEQ(Z + Z^2). """

        spec = Specification()
        z, C = Variable(), Variable()
        spec.add(C, Seq(z + z**2))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.618034527351341, 5) # golden ratio

    def test_singular_partitions(self):
        """ Singular tuning of partitions
            P = MSET(SEQ_{k >= 1}(Z))."""

        params = Params(Type.ALGEBRAIC)
        params.max_iters = 1000 # required

        spec = Specification()
        z, P = Variable(), Variable()
        spec.add(P, MSet(z * Seq(z)))

        spec.run_singular_tuner(z, params)
        self.assertAlmostEqual(z.value, 0.999992520391430, 5)

    def test_minus_constant(self):

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, Seq(2*z) - 1)

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.5, 5)

    def test_minus_constant2(self):

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, z - 2 * T)

        try:
            spec.run_singular_tuner(z)
        except ValueError:
            self.assertTrue(z.value is None)

    def test_binary_necklaces(self):
        """ Singular tuning of neckleces build using two kinds of beads.
            N = CYC(Z + Z)."""

        spec = Specification()
        z, N = Variable(), Variable()
        spec.add(N, Cyc(z + z))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.5, 5)

    def test_cyclic_compositions(self):
        """ Singular tuning of cyclic compositions.
            C = CYC(Z * SEQ(Z))."""

        spec = Specification()
        z, C = Variable(), Variable()
        spec.add(C, Cyc(z * Seq(z)))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.5, 5)

    def test_unlabelled_functional_graphs(self):
        """ Singular tuning of unlabelled functional graphs.
            F = MSet(K)
            K = CYC(U)
            U = Z * MSet(U)."""

        spec = Specification()
        z, F = Variable(), Variable()
        K, U = Variable(), Variable()

        spec.add(F, MSet(K))
        spec.add(K, Cyc(U))
        spec.add(U, z * MSet(U))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.3383218568992077, 5)

    def test_ternary_trees(self):
        """ Singular ternary trees.
            T = 1 + Z * Seq_{= 3}(Z)."""

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, 1 + z * Seq(T, eq(3)))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.148148148148149, 5)

class MeanTuner(unittest.TestCase):

    def test_motzkin_trees(self):
        """ Tuning of Motzkin trees
            M = Z + Z * M + Z * M ** 2.
            (expected size around 1000)"""

        spec = Specification()
        z, M = Variable(1000), Variable()
        spec.add(M, z + z * M + z * M ** 2)

        params = Params(Type.ALGEBRAIC)
        spec.run_tuner(M, params)

        self.assertAlmostEqual(z.value, 0.333333083333287)
        self.assertAlmostEqual(M.value, 0.998501123876053)

    def test_lambda_terms(self):
        """ Tuning of lambda terms
            L = D + Z * L + Z * L ** 2
            D = Z + Z * D"""

        spec = Specification()
        z = Variable(10000) # size
        v = Variable(3120)  # variables
        u = Variable(312)   # successors

        L, D = Variable(), Variable()

        spec.add(L, D + z * L + z * L ** 2)
        spec.add(D, v * z + u * z * D)

        params = Params(Type.ALGEBRAIC)
        spec.run_tuner(L, params)

        self.assertAlmostEqual(z.value, 0.356007431874485)
        self.assertAlmostEqual(L.value, 0.904320092780514)

    def test_lambda_terms2(self):
        """ Tuning of lambda terms
            L = D + Z * L + Z * L ** 2
            D = Z + Z * D"""

        spec = Specification()
        z = Variable()     # size
        u = Variable(0.4)  # abstractions

        L, D = Variable(), Variable()

        spec.add(L, D + u * z * L + z * L ** 2)
        spec.add(D, z + z * D)

        params = Params(Type.ALGEBRAIC)
        spec.run_singular_tuner(z, params)

        self.assertAlmostEqual(z.value, 0.244827373512259)
        self.assertAlmostEqual(u.value, 1.78303233505684)
        self.assertAlmostEqual(L.value, 1.15073912781323)

if __name__ == '__main__':
    unittest.main()
