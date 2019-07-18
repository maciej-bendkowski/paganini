import unittest

from paganini.utils import *
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

        self.assertAlmostEqual(z.value, 0.338322112871298, 5)
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

    def test_ternary_trees(self):
        """ Singular ternary trees.
            T = 1 + Z * Seq_{= 3}(Z)."""

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, 1 + z * Seq(T, eq(3)))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.148148148148149, 5)

    def test_otter_trees(self):
        """ Singular Otter trees.
            T = 1 + Z * MSet_{ = 2}(T)."""

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, 1 + z * MSet(T, eq(2)))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.4026975, 5)

    def test_otter_trees2(self):
        """ Singular Otter trees.
            T = 1 + Z * MSet_{ = 3}(T)."""

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, 1 + z * MSet(T, eq(3)))

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.355181762886292, 5)

    def test_custom_singular_btrees(self):
        """ Singular, custom btrees."""

        spec = Specification()
        z, a, b, T = Variable(), Variable(0.5), Variable(0.5), Variable()
        spec.add(T, z * (a + b) + T * T)

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.249999999878295, 5)


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

    def test_cyclic_compositions2(self):
        """ Tuning of bounded cyclic compositions.
            C = CYC_{= 12}(Z * SEQ(Z))."""

        spec = Specification()
        z, C = Variable(20), Variable()
        spec.add(C, UCyc(z * Seq(z), eq(12)))

        spec.run_tuner(C)
        self.assertAlmostEqual(z.value, 0.405765659263783, 5)

    def test_urns(self):
        """ Tuning of urns.
            U = Set(Z)."""

        spec = Specification()
        z, U = Variable(18), Variable()
        spec.add(U, Set(2 * z))

        spec.run_tuner(U)
        self.assertAlmostEqual(z.value, 9, 5)

    def test_seq_urns2(self):
        """ Tuning of sequences of bounded urns.
            U = Seq(Set_{= 3}(Z))."""

        spec = Specification()
        z, U = Variable(5), Variable()
        spec.add(U, Seq(Set(z, eq(3))))

        spec.run_tuner(U)
        self.assertAlmostEqual(z.value, 1.55361625297693, 5)

    def test_circular_graphs(self):
        """ Tuning of circular graphs.
            C = Cyc(z)."""

        spec = Specification()
        z, C = Variable(10), Variable()
        spec.add(C, Cyc(z))

        spec.run_tuner(C)
        self.assertAlmostEqual(z.value, 1.12975951282490, 5)

    def test_alignments(self):
        """ Tuning of alignments.
            O = Seq(Cyc(z))."""

        spec = Specification()
        z, O = Variable(10), Variable()
        spec.add(O, Seq(Cyc(z)))

        spec.run_tuner(O)
        self.assertAlmostEqual(z.value, 0.578097783364826, 5)

    def test_permutations(self):
        """ Tuning of permutations.
            P = Set(Cyc(z))."""

        spec = Specification()
        z, P = Variable(666), Variable()
        spec.add(P, Set(Cyc(z)))

        spec.run_tuner(P)
        self.assertAlmostEqual(z.value, 1.28397586928450, 5)

    def test_set_permutations(self):
        """ Tuning of set permutations.
            P = Set(Set_{>= 1}(z))."""

        spec = Specification()
        z, P = Variable(32), Variable()

        # approximation Set_{>= 1} = sum_{k = 1}^K Set_{= k}.
        s = sum([Set(z, eq(k)) for k in range(1, 20)])
        spec.add(P, Set(s))

        spec.run_tuner(P)
        self.assertAlmostEqual(z.value, 5.25734205219187, 5)

    def test_set_surjections(self):
        """ Tuning of set surjections.
            S = Seq(Set_{>= 1}(z))."""

        spec = Specification()
        z, S = Variable(32), Variable()

        # approximation Set_{>= 1} = sum_{k = 1}^K Set_{= k}.
        s = sum([Set(z, eq(k)) for k in range(1, 20)])
        spec.add(S, Seq(s))

        spec.run_tuner(S)
        self.assertAlmostEqual(z.value, 0.672353796989521, 5)

    def test_arrangements(self):
        """ Tuning of arrangements.
            A = U P
            U = Set(z)
            P = Seq(z)."""

        spec = Specification()
        z, A, = Variable(1024), Variable()
        U, P  = Variable(), Variable()

        spec.add(A, U * P)
        spec.add(U, Set(z))
        spec.add(P, Seq(z))

        spec.run_tuner(A)
        self.assertAlmostEqual(z.value, 0.999023438431325, 5)

    def test_derangements(self):
        """ Tuning of derangements.
            D = Set(Cyc_{> 3}(z))."""

        spec = Specification()
        z, D = Variable(10), Variable()

        # approximation Cyc_{> 1} = sum_{k = 2}^K Cyc_{= k}.
        cs = sum([Cyc(z, eq(k)) for k in range(4, 24)])

        spec.add(D, Set(cs))

        spec.run_tuner(D)
        self.assertAlmostEqual(z.value, 1.18802573842469, 5)

    def test_cayley_trees(self):
        """ Tuning of Cayley trees.
            T = Z Set(T)."""

        spec = Specification()
        z, T = Variable(1024), Variable()
        spec.add(T, z * Set(T))

        spec.run_tuner(T)
        self.assertAlmostEqual(z.value, 0.367879265638609, 5)

class UtilsTuner(unittest.TestCase):

    def test_partition_sequences(self):
        """ Checks that each of the generated partition-sequences
        has proper length and structure (sums up to its length)."""

        for n in range(2, 20):
            for ps in partition_sequences(n):
                self.assertEqual(len(ps), n)

                total = 0
                for i, k in enumerate(ps):
                    total += (i + 1) * k

                self.assertEqual(total, n)

class ExpressionsTest(unittest.TestCase):

    def test_related_expressions(self):
        x, y, z = Variable(), Variable(), Variable()

        self.assertTrue(x.related(x))
        self.assertFalse(x.related(y))
        self.assertFalse(x.related(x * x))
        self.assertTrue((x * y).related(y * x))
        self.assertFalse((x * y * z).related(y * x))

if __name__ == '__main__':
    unittest.main()
