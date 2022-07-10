import unittest

from paganini.expressions import *
from paganini.specification import *
from paganini.utils import *


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
        self.assertGreater(len(spec.ddg(B)), 0)

    def test_singular_motzkin_trees(self):
        """ Singular tuning of Motzkin trees
            M = Z * SEQ<=2(M). """

        spec = Specification()
        z, M = Variable(), Variable()

        seqM = Seq(M, leq(2))
        spec.add(M, z * seqM)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.333333333333334)
        self.assertAlmostEqual(M.value, 1.0)

        self.assertEqual(len(spec.ddg(M)), 0)
        self.assertGreater(len(spec.ddg(seqM)), 0)

    def test_singular_motzkin_trees2(self):
        """ Singular tuning of Motzkin trees
            M = Z + Z * M + Z * M^2. """

        spec = Specification()
        z, M = Variable(), Variable()
        spec.add(M, z + z * M + z * M ** 2)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.333333333333334)
        self.assertAlmostEqual(M.value, 1.0)

        self.assertGreater(len(spec.ddg(M)), 0)

    def test_singular_trees(self):
        """ Singular tuning of plane trees
            T = Z * SEQ(T)."""

        spec = Specification()
        z, T = Variable(), Variable()

        seqT = Seq(T)
        spec.add(T, z * seqT)

        spec.run_singular_tuner(z)

        self.assertAlmostEqual(z.value, 0.25)
        self.assertAlmostEqual(T.value, 0.5)

        self.assertEqual(len(spec.ddg(T)), 0)
        self.assertGreater(len(spec.ddg(seqT)), 0)

    def test_singular_lambda_terms(self):
        """ Singular tuning of plain lambda terms
            L = Z * SEQ(Z) + Z * L + Z * L^2."""

        spec = Specification()
        z, L, D = Variable(), Variable(), Variable()
        spec.add(L, D + z * L + z * L ** 2)
        spec.add(D, z + z * D)

        spec.run_singular_tuner(z, method=Method.FORCE)

        self.assertAlmostEqual(z.value, 0.295597742522085)
        self.assertAlmostEqual(L.value, 1.19148788395312)

        self.assertGreater(len(spec.ddg(L)), 0)
        self.assertGreater(len(spec.ddg(D)), 0)

    def test_singular_lambda_terms2(self):
        """ Singular tuning of plain lambda terms
            L = Z * SEQ(Z) + Z * L + Z * L^2."""

        spec = Specification()
        z, L = Variable(), Variable()
        spec.add(L, z * Seq(z) + z * L + z * L ** 2)

        spec.run_singular_tuner(z, method=Method.FORCE)

        self.assertAlmostEqual(z.value, 0.295597742522085)
        self.assertAlmostEqual(L.value, 1.19148788395312)

        self.assertGreater(len(spec.ddg(L)), 0)

    def test_singular_polya_trees(self):
        """ Singular tuning of Polya trees
            T = Z * MSet(T)."""

        spec = Specification()
        z, T = Variable(200000), Variable()

        msetT = MSet(T)
        spec.add(T, z * msetT)

        spec.run_tuner(T)

        self.assertAlmostEqual(z.value, 0.338322112871298, 5)
        self.assertAlmostEqual(T.value, 0.999993919977244)

        self.assertEqual(len(spec.ddg(T)), 0)
        self.assertGreater(len(spec.ddg(msetT)), 0)

        diags = spec.diagonals()
        self.assertGreater(len(diags[T]), 0)
        self.assertAlmostEqual(diags[T][1].value, T.value)

    def test_singular_custom_trees(self):
        """ Singular tuning of some custom trees defined by
            T = Z + Z * SEQ_>=2(T)."""

        params = Params(Type.ALGEBRAIC)
        params.max_iters = 100  # required

        spec = Specification()
        z, T = Variable(), Variable()
        seqT = Seq(T, geq(2))
        spec.add(T, z + z * seqT)

        spec.run_singular_tuner(z, params)

        self.assertAlmostEqual(z.value, 0.333333333333335)
        self.assertAlmostEqual(T.value, 0.499999999999993)

        # note: = seq(T)_{>=2} = T^2 Seq(T)
        self.assertGreater(len(spec.ddg(T)), 0)
        self.assertEqual(len(spec.ddg(seqT)), 0)

    def test_binary_words(self):
        """ Singular tuning of binary words.
            B = SEQ(Z + Z). """

        spec = Specification()
        z, B = Variable(), Variable()
        seqZ = Seq(z + z)
        spec.add(B, seqZ)

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.5, 5)

        self.assertEqual(len(spec.ddg(seqZ)), 0)  # note: trivial distribution.
        self.assertEqual(len(spec.ddg(B)), 0)

    def test_compositions(self):
        """ Singular tuning of all compositions.
            C = SEQ(Z * SEQ(Z)). """

        params = Params(Type.RATIONAL)
        params.max_iters = 10000  # required

        spec = Specification()
        z, C = Variable(), Variable()
        seq = Seq(z * Seq(z))
        spec.add(C, seq)

        spec.run_singular_tuner(z, params, Method.FORCE)
        self.assertAlmostEqual(z.value, 0.5, 3)

        self.assertEqual(len(spec.ddg(C)), 0)
        self.assertEqual(len(spec.ddg(seq)), 0)

    def test_compositions_with_restricted_summands(self):
        """ Singular tuning of compositions with restricted summands in {1,2}.
            C = SEQ(Z + Z^2). """

        spec = Specification()
        z, C = Variable(500000), Variable()
        spec.add(C, Seq(z + z ** 2))

        spec.run_tuner(C)
        self.assertAlmostEqual(z.value, 0.618034527351341, 5)  # golden ratio

    def test_singular_partitions(self):
        """ Singular tuning of partitions
            P = MSET(SEQ_{k >= 1}(Z))."""

        params = Params(Type.ALGEBRAIC)

        spec = Specification()
        z, P = Variable(200000), Variable()
        mset = MSet(z * Seq(z))
        spec.add(P, mset)

        spec.run_tuner(P, params)
        self.assertAlmostEqual(z.value, 0.997247023094167, 5)

        self.assertEqual(len(spec.ddg(P)), 0)
        self.assertEqual(len(spec.ddg(mset)), 0)  # note: trivial distribution

    def test_minus_constant(self):

        spec = Specification()
        z, T = Variable(50000), Variable()
        spec.add(T, Seq(2 * z) - 1)

        spec.run_tuner(T)
        self.assertAlmostEqual(z.value, 0.5, 4)

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

        self.assertGreater(len(spec.ddg(T)), 0)

    def test_otter_trees(self):
        """ Singular Otter trees.
            T = 1 + Z * MSet_{ = 2}(T)."""

        with self.assertRaises(ValueError):
            spec = Specification()
            z, T = Variable(), Variable()
            spec.add(T, 1 + z * MSet(T, eq(2)))
            spec.run_singular_tuner(z)

        spec = Specification()
        z, T = Variable(1000000), Variable()
        spec.add(T, 1 + z * MSet(T, eq(2)))

        spec.run_tuner(T)
        self.assertAlmostEqual(z.value, 0.4026975, 5)

        self.assertGreater(len(spec.ddg(T)), 0)

    def test_otter_trees2(self):
        """ Singular Otter trees.
            T = 1 + Z * MSet_{ = 3}(T)."""

        with self.assertRaises(ValueError):
            spec = Specification()
            z, T = Variable(), Variable()
            spec.add(T, 1 + z * MSet(T, eq(3)))
            spec.run_singular_tuner(z)

        spec = Specification()
        z, T = Variable(), Variable()
        spec.add(T, 1 + z * MSet(T, eq(3)))

        # note: Method.FORCE is not necessary here.
        spec.run_singular_tuner(z, method=Method.FORCE)
        self.assertAlmostEqual(z.value, 0.355181762886292, 5)

        self.assertGreater(len(spec.ddg(T)), 0)

    def test_custom_singular_btrees(self):
        """ Singular, custom btrees."""

        spec = Specification()
        z, a, b, T = Variable(), Variable(0.5), Variable(0.5), Variable()
        spec.add(T, z * (a + b) + T * T)

        spec.run_singular_tuner(z)
        self.assertAlmostEqual(z.value, 0.249999999878295, 5)

        self.assertGreater(len(spec.ddg(T)), 0)


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

        self.assertGreater(len(spec.ddg(M)), 0)

    def test_lambda_terms(self):
        """ Tuning of lambda terms
            L = D + Z * L + Z * L ** 2
            D = Z + Z * D"""

        spec = Specification()
        z = Variable(10000)  # size
        v = Variable(3120)  # variables
        u = Variable(312)  # successors

        L, D = Variable(), Variable()

        spec.add(L, D + z * L + z * L ** 2)
        spec.add(D, v * z + u * z * D)

        params = Params(Type.ALGEBRAIC)
        spec.run_tuner(L, params)

        self.assertAlmostEqual(z.value, 0.356007431874485)
        self.assertAlmostEqual(L.value, 0.904320092780514)

        self.assertGreater(len(spec.ddg(L)), 0)
        self.assertGreater(len(spec.ddg(D)), 0)

    def test_lambda_terms2(self):
        """ Tuning of lambda terms
            L = D + Z * L + Z * L ** 2
            D = Z + Z * D"""

        spec = Specification()
        z = Variable(100000)  # size
        u = Variable(40000)  # abstractions

        L, D = Variable(), Variable()

        spec.add(L, D + u * z * L + z * L ** 2)
        spec.add(D, z + z * D)

        params = Params(Type.ALGEBRAIC)
        spec.run_tuner(L, params)

        self.assertAlmostEqual(z.value, 0.244827141130008, 5)
        self.assertAlmostEqual(u.value, 1.7830323350568, 5)
        self.assertAlmostEqual(L.value, 1.1507391278132, 4)

        self.assertGreater(len(spec.ddg(L)), 0)
        self.assertGreater(len(spec.ddg(D)), 0)

    def test_cyclic_compositions2(self):
        """ Tuning of bounded cyclic compositions.
            C = CYC_{= 12}(Z * SEQ(Z))."""

        spec = Specification()
        z, C = Variable(20), Variable()
        ucyc = UCyc(z * Seq(z), eq(12))
        spec.add(C, ucyc)

        spec.run_tuner(C)
        self.assertAlmostEqual(z.value, 0.405765659263783, 5)

        self.assertGreater(len(spec.ddg(ucyc)), 0)

    def test_urns(self):
        """ Tuning of urns.
            U = Set(Z)."""

        spec = Specification()
        z, U = Variable(18), Variable()
        setz = Set(2 * z)
        spec.add(U, setz)

        spec.run_tuner(U)
        self.assertAlmostEqual(z.value, 9, 5)

        self.assertEqual(len(spec.ddg(setz)), 0)  # note: trivial distribution

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
        cycz = Cyc(z)
        spec.add(C, cycz)

        spec.run_tuner(C)
        self.assertAlmostEqual(z.value, 1.12975951282490, 5)

        self.assertGreater(len(spec.ddg(cycz)), 0)

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
        U, P = Variable(), Variable()

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

    def test_near_singular_forests_broken(self):
        """ Singular forests with 'unreachable' trees. """

        spec = Specification()
        z, Tree, Forest = Variable(1000000), Variable(), Variable()

        spec.add(Forest, Seq(Tree))
        spec.add(Tree, z + Tree * Tree)

        spec.run_tuner(Forest)

        self.assertAlmostEqual(z.value, 0.25, 5)
        self.assertAlmostEqual(Tree.value, 0.5, 5)
        self.assertAlmostEqual(Forest.value, 2, 5)

    def test_near_singular_forests(self):
        """ Singular forests with 'unreachable' trees. """

        spec = Specification()
        z, Tree, Forest = Variable(), Variable(), Variable()

        spec.add(Forest, Seq(Tree))
        spec.add(Tree, z + Tree * Tree)

        with self.assertRaises(ValueError):
            spec.run_singular_tuner(z)

    def test_finite_binary_words(self):
        """ Tuning of binary words.
            B = SEQ(Z + Z). """

        spec = Specification()
        z, B = Variable(100000), Variable()
        spec.add(B, Seq(z + z))

        spec.run_tuner(B)
        self.assertAlmostEqual(z.value, 0.5, 5)

    def test_non_reachable_states(self):
        """ Specification with unreachable states. """
        spec = Specification()
        z, T1, T2 = Variable(5000), Variable(), Variable()

        spec.add(T1, 1 + z * T1 ** 2)
        spec.add(T2, 1 + z * T1 ** 2 + z * T2 ** 3)

        with self.assertRaises(ValueError):
            spec.run_tuner(T1)

    def test_non_reachable_states2(self):
        """ Specification with unreachable states.
        Tuned in the 'proper' way. """

        spec = Specification()
        z, T1, T2 = Variable(100000), Variable(), Variable()

        spec.add(T1, 1 + z * T1 ** 2)
        spec.add(T2, 1 + z * T1 ** 2 + z * T2 ** 3)

        spec._check_singular_tuner()
        problem = spec.run_tuner(T2)
        self.assertAlmostEqual(z.value, 0.112382934442027, 5)

    def test_not_well_founded(self):
        spec = Specification()
        z, B = Variable(2000), Variable(100)

        spec.add(B, 1 + B ** 2)
        problem = spec.run_tuner(B)
        self.assertEqual(problem, float("inf"))


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
    def assertExprEqual(self, x, y):
        self.assertIsInstance(x, Expr)
        self.assertIsInstance(y, Expr)
        self.assertEqual(x.coeff, y.coeff)
        self.assertEqual(x.variables, y.variables)

    def test_related_expressions(self):
        x, y, z = Variable(), Variable(), Variable()

        self.assertTrue(x.related(x))
        self.assertFalse(x.related(y))
        self.assertFalse(x.related(x * x))
        self.assertTrue((x * y).related(y * x))
        self.assertFalse((x * y * z).related(y * x))

    def test_exponentiation(self):
        x = Variable()

        self.assertExprEqual(x ** 0, Expr(1))
        self.assertExprEqual(x ** 1, x)
        self.assertExprEqual(x ** 3, x * x * x)

    def test_polynomial_variable_check(self):

        x = Variable()
        px = Polynomial.cast(x)
        self.assertTrue(px.is_variable())

        py = Polynomial.cast(x + 2 * x)
        self.assertFalse(py.is_variable())


if __name__ == "__main__":
    unittest.main()
