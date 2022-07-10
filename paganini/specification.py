from enum import Enum
from fractions import Fraction
from math import factorial

import cvxpy
import networkx as nx
import numpy as np
import sympy
from optas.construct import construct_sample_ky_encoding
from optas.divergences import KERNELS
from optas.opt import get_optimal_probabilities

from paganini.expressions import *
from paganini.utils import partition_sequences, phi

# define the namespace to avoid polluting with foreign packages
__all__ = (
    "Seq",
    "UCyc",
    "MSet",
    "Set",
    "Cyc",
    "Operator",
    "Constraint",
    "Type",
    "Params",
    "Method",
    "Specification",
    "leq",
    "geq",
    "eq",
)


class Seq(Variable):
    """ Sequence variables."""

    def __init__(self, expression, constraint=None):
        super(Seq, self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE  # make sure its a type variable
        self.constraint = Constraint.normalise(constraint)

    def register(self, spec):
        """ Unfolds the Seq definition and registers it in the given system."""

        if self.constraint.operator == Operator.UNBOUNDED:
            # note: Seq(expr) = 1 + expr * Seq(expr).
            spec.add(self, 1 + self.inner_expressions * self)

        elif self.constraint.operator == Operator.LEQ:
            # note: Seq(expr)_{<= k} = 1 + expr + expr^2 + ... + expr^k.
            expressions = Expr(1)
            for k in range(1, self.constraint.value + 1):
                expressions = expressions + self.inner_expressions ** k

            spec.add(self, expressions)

        elif self.constraint.operator == Operator.EQ:
            # note: Seq(expr)_{= k} = expr ** k
            spec.add(self, self.inner_expressions ** self.constraint.value)

        else:  # constraint.operator == Operator.GEQ
            # note: Seq(expr)_{>= k} = expr^k + expr^{k+1} + ...
            #                        = expr^k (1 + expr^2 + expr^3 + ...)
            #                        = expr^k Seq(expr).

            seq = Seq(self.inner_expressions)
            spec.add(seq, 1 + self.inner_expressions * seq)
            spec.add(self, self.inner_expressions ** self.constraint.value * seq)


class MSet(Variable):
    """ MSet variables."""

    def __init__(self, expression, constraint=None):
        super(MSet, self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE  # make sure its a type variable
        self.constraint = Constraint.normalise(constraint)

        if (
            self.constraint.operator == Operator.LEQ
            or self.constraint.operator == Operator.GEQ
        ):
            raise AttributeError("Unsupported constraint.")

    def register(self, spec):
        """ Unfolds the MSet definition and registers it in the given system."""
        spec._diagonal_variable(self, 1)


class UCyc(Variable):
    """ Unlabelled Cyc variables."""

    def __init__(self, expression, constraint=None):
        super(UCyc, self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE  # make sure its a type variable
        self.constraint = Constraint.normalise(constraint)

        if self.constraint.operator != Operator.EQ:
            raise AttributeError("Unsupported constraint.")

    def register(self, spec):
        """ Unfolds the UCyc definition and registers it in the given system."""
        spec._diagonal_variable(self, 1)


class Set(Variable):
    """ Labelled Set variables."""

    def __init__(self, expression, constraint=None):
        super(Set, self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE  # make sure its a type variable
        self.constraint = Constraint.normalise(constraint)

        if (
            self.constraint.operator == Operator.LEQ
            or self.constraint.operator == Operator.GEQ
        ):
            raise AttributeError("Unsupported constraint.")

    def register(self, spec):
        """ Unfolds the Set definition and registers it in the given system."""
        spec._diagonal_variable(self, 1)


class Cyc(Variable):
    """ Labelled Cyc variables."""

    def __init__(self, expression, constraint=None):
        super(Cyc, self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE  # make sure its a type variable
        self.constraint = Constraint.normalise(constraint)

        if (
            self.constraint.operator == Operator.LEQ
            or self.constraint.operator == Operator.GEQ
        ):
            raise AttributeError("Unsupported constraint.")

    def register(self, spec):
        """ Unfolds the Cyc definition and registers it in the given system."""
        spec._diagonal_variable(self, 1)


class Operator(Enum):
    """ Enumeration of supported constraint signs."""

    LEQ = 1  # less or equal
    EQ = 2  # equal
    GEQ = 3  # greater or equal
    UNBOUNDED = 4  # unbounded operator


class Constraint:
    """ Supported constraints for classes such as SEQ."""

    def __init__(self, operator, value):
        self.operator = operator
        self.value = value

    @staticmethod
    def normalise(constraint=None):

        if constraint is None:
            return Constraint(Operator.UNBOUNDED, 0)
        else:
            return constraint


def leq(n):
    """ Creates a less or equal constraint for the given input."""
    assert n >= 0, "Negative constraints are not supported."
    return Constraint(Operator.LEQ, n)


def eq(n):
    """ Creates an equal constraint for the given input."""
    assert n > 0, "Non-positive constraints are not supported."
    return Constraint(Operator.EQ, n)


def geq(n):
    """ Creates a greater or equal constraint for the given input."""
    assert n >= 0, "Negative constraints are not supported."
    return Constraint(Operator.GEQ, n)


class Type(Enum):
    """ Enumeration of supported system types."""

    ALGEBRAIC = 1
    RATIONAL = 2


class Params:
    """ CVXPY solver parameters initalised with some defaults."""

    def __init__(self, sys_type):
        self.verbose = False
        if sys_type == Type.RATIONAL:
            self.sys_type = Type.RATIONAL
            self.solver = cvxpy.SCS
            self.max_iters = 2500
            self.eps = 1.0e-20
            self.norm = 40
        else:
            self.sys_type = Type.ALGEBRAIC
            self.solver = cvxpy.ECOS
            self.max_iters = 250
            self.feastol = 1.0e-20
            self.abstol = 1.0e-20
            self.reltol = 1.0e-20


class Method(Enum):
    """ Enumeration of supported method flags. See specification.run_singular_tuner."""

    STRICT = 1
    FORCE = 2  # skip theoretical checks.


class Specification:
    """ Symbolic system specifications."""

    def __init__(self, series_truncate=20):
        """ Creates a new specification. The optional `series_truncate`
        parameter controls the truncation threshold for infinite series,
        intrinsic to some more involved constructions, such as multisets or
        cycles."""

        self._index_counter = 0
        self._equations = {}
        self._graph = nx.DiGraph()

        self._all_variables = {}
        self._tuned_variables = set()

        self._seq_variables = set()
        self._mset_variables = set()
        self._ucyc_variables = set()

        self._set_variables = set()
        self._cyc_variables = set()

        # diagonals.
        self._diag = {}
        self._msets = {}
        self._sets = {}

        self._series_truncate = series_truncate

    def diagonals(self):
        """ Map from Pólya variables to their "diagonals".
        For instance,

        >>> T, z = Variable(), Variable()
        >>> spec = Specification()
        >>> spec.add(T, z * MSet(T))
        >>> spec.run_singular_tuner(z, method=Method.FORCE)

        >>> for j, var in spec.diagonals()[T].items():
        >>>     print(f"T(z^{j}) = {var.value}")
        """
        return self._diag

    def __repr__(self):
        return "\n".join(
            [
                v.__repr__() + " = " + self._equations[v].__repr__()
                for v in self._equations
            ]
        )

    def _next_idx(self):
        n = self._index_counter
        self._index_counter += 1
        return n

    @property
    def discharged_variables(self):
        """ Number of variables discharged in the system."""
        return self._index_counter

    def _register_variable(self, v):
        """ Registers the given variable in the specification."""

        if v.idx is not None:
            return  # nothing to do.

        v.idx = self._next_idx()
        self._all_variables[v.idx] = v

        if v.tuning_param is not None:
            self._tuned_variables.add(v)

        # case over the variable's type
        if isinstance(v, Seq):
            self._seq_variables.add(v)
            self._register_expressions(v.inner_expressions)

        elif isinstance(v, MSet):
            self._mset_variables.add(v)
            self._register_expressions(v.inner_expressions)

        elif isinstance(v, UCyc):
            self._ucyc_variables.add(v)
            self._register_expressions(v.inner_expressions)

        elif isinstance(v, Set):
            self._set_variables.add(v)
            self._register_expressions(v.inner_expressions)

        elif isinstance(v, Cyc):
            self._cyc_variables.add(v)
            self._register_expressions(v.inner_expressions)

    def _register_expression(self, expression):
        assert isinstance(expression, Expr), "Expected expression."
        for v in expression.variables:
            self._register_variable(v)

    def _register_expressions(self, expressions):
        if isinstance(expressions, Expr):
            self._register_expression(expressions)
        else:
            for expr in expressions:
                self._register_expression(expr)

    def _type_variable(self):
        """ Discharges a fresh variable."""
        v = Variable()
        v.type = VariableType.TYPE
        self._register_variable(v)
        return v

    def _build_dg_variable(self, source, var):
        assert isinstance(var, Variable), "Expected variable."
        self._add_dependency(source, var)
        if hasattr(var, "inner_expressions"):
            self._build_dg(var, var.inner_expressions)

        # case over the variable's type adding dependencies
        # between 'var' and its (truncated) series variables.
        if isinstance(var, MSet) and var in self._msets:
            self._build_dg(var, Polynomial.sum(self._msets[var]))

        if isinstance(var, Set) and var in self._sets:
            self._build_dg(var, Polynomial.sum(self._sets[var]))

    def _build_dg_expr(self, source, expr):
        assert isinstance(expr, Expr), "Expected expression."
        for v in expr.variables:
            self._build_dg_variable(source, v)

    def _build_dg(self, source, expressions):
        if isinstance(expressions, Expr):
            self._build_dg_expr(source, expressions)
        else:
            for expr in expressions:
                self._build_dg_expr(source, expr)

    def _build_dependency_graph(self):

        if len(self._graph) > 0:
            return  # do nothing if already built

        for v in self._equations:
            rhs = self._equations[v]
            self._build_dg(v, rhs)

    def check_type(self):
        """ Checks if the system is algebraic or rational.
        Note: the current method is heuristic."""

        if (
            len(self._seq_variables) > 0
            or len(self._mset_variables) > 0
            or len(self._ucyc_variables) > 0
            or len(self._cyc_variables) > 0
            or len(self._set_variables) > 0
        ):
            return Type.ALGEBRAIC

        for expressions in self._equations.values():
            for exp in expressions:
                if isinstance(exp, Expr):
                    for v, e in exp.variables.items():
                        if v.is_type_variable and e > 1:
                            return Type.ALGEBRAIC

        return Type.RATIONAL

    def _init_params(self, params=None):
        if params is None:
            # some defaults
            sys_type = self.check_type()
            return Params(sys_type)
        else:
            return params

    def _diagonal_expr(self, expr, d=1):
        """ Extends self._diag_variable to expressions."""
        assert isinstance(expr, Expr), "Expression required."

        variables = {}
        for v in expr.variables:
            if v.is_type_variable:
                # substitute the power variable.
                x = self._diagonal_variable(v, d)
                variables[x] = expr.variables[v]
            else:
                # increase the exponent.
                variables[v] = d * expr.variables[v]

        return Expr(expr.coeff, variables)

    def _diagonal_variable(self, var, d=1):
        """ Given :math:`T(z)` creates its kth diagonal, i.e. :math:`T(z^k)`."""

        assert d > 0, "Non-positive diagonal parameter."
        assert var.is_type_variable, "Requested diagonal of non-type variable."

        if var not in self._diag:
            self._diag[var] = {}  # degree -> expression

        # check the variable cache.
        if d in self._diag[var]:
            return self._diag[var][d]

        # trivial diagonal case.
        if d == 1 and var in self._equations:
            self._diag[var][d] = var
            return var

        var_d = self._type_variable() if d > 1 else var
        self._diag[var][d] = var_d  # memorise var[d]

        if var in self._equations:
            monomials = self._equations[var]
            self.add(var_d, [self._diagonal_expr(e, d) for e in monomials])
            return var_d

        if var in self._mset_variables:

            if var.constraint.operator == Operator.EQ:
                # note: constrained MSet_{ = k}.

                series, k = [], var.constraint.value
                for n in partition_sequences(k):
                    product = Polynomial(Expr(1))
                    for i in range(1, k + 1):  # note the truncation.
                        if i * d <= self._series_truncate:
                            c = 1 / (factorial(n[i - 1]) * (i ** n[i - 1]))
                            expr = Polynomial(
                                [
                                    self._diagonal_expr(e, i * d)
                                    for e in var.inner_expressions
                                ]
                            )

                            product *= c * expr ** n[i - 1]

                    if not product.is_one():
                        series.append(product)

                self.add(var_d, Polynomial.sum(series))
                return var_d

            else:
                self._msets[var_d] = []
                for k in range(d, self._series_truncate + 1, d):
                    self._msets[var_d].append(
                        Polynomial(
                            [self._diagonal_expr(e, k) for e in var.inner_expressions]
                        )
                    )

                return var_d

        elif var in self._ucyc_variables:
            if var.constraint.operator == Operator.EQ:
                # note: constrained Cyc_{ = k}.

                series, k = [], var.constraint.value
                for i in range(1, k + 1):
                    if k % i == 0:
                        expr = Polynomial(
                            [
                                self._diagonal_expr(e, i * d)
                                for e in var.inner_expressions
                            ]
                        )

                        series.append((phi(i) / k) * expr ** (k // i))

                self.add(var_d, Polynomial.sum(series))
                return var_d

        elif var in self._set_variables:

            if var.constraint.operator == Operator.EQ:
                # note: constrained Set_{ = k}.

                k = var.constraint.value
                expr = Polynomial(
                    [self._diagonal_expr(e, d) for e in var.inner_expressions]
                )

                self.add(var_d, (1 / (factorial(k))) * expr ** k)
                return var_d

            else:
                self._sets[var_d] = [
                    Polynomial(
                        [self._diagonal_expr(e, d) for e in var.inner_expressions]
                    )
                ]

                return var_d

        elif var in self._cyc_variables:

            if var.constraint.operator == Operator.EQ:
                # note: constrained Cyc_{ = k}.

                k = var.constraint.value
                expr = Polynomial(
                    [self._diagonal_expr(e, d) for e in var.inner_expressions]
                )

                self.add(var_d, (1 / k) * expr ** k)
                return var_d

            else:
                series = []
                for k in range(1, self._series_truncate + 1):
                    series.append(
                        1
                        / k
                        * Polynomial(
                            [self._diagonal_expr(e, d) for e in var.inner_expressions]
                        )
                        ** k
                    )

                self.add(var_d, Polynomial.sum(series))
                return var_d

    def _add_dependency(self, source, target):
        """ Connects variable 'source' to its dependency 'target'. If source is
        None, no dependency is formed. Note: We assume that both variables are
        already registered in the specification.
        """

        if source is not None:
            self._graph.add_node(source.idx)  # no multiple nodes.
            self._graph.add_node(target.idx)
            self._graph.add_edge(source.idx, target.idx)  # no multiple edges.

    def add(self, var, expression):
        """ Includes the given definition in the specification."""
        expression = Polynomial.cast(expression)
        var.type = VariableType.TYPE  # make var a type variable
        self._equations[var] = expression

        # register variables in the system.
        self._register_variable(var)
        self._register_expressions(expression)

    def _is_synonym_variable(self, v):

        if v not in self._equations:
            return False

        if v.idx not in self._graph:
            return False

        if self._graph.in_degree(v.idx) != 0:
            return False

        return self._equations[v].is_variable()

    def _compose_constraints(self, variables, n):
        """ Composes optimisation constraints."""

        constraints = []
        for v in self._equations:
            rhs = self._equations[v]  # right-hand side.
            matrix, coeffs, constant_term = rhs.specification(n)

            exponents = matrix @ variables + coeffs
            if self._is_synonym_variable(v):

                x = self._equations[v]._expressions[0]
                constraints.append(variables[x.idx] == variables[v.idx])
            else:
                constraints.append(
                    variables[v.idx] >= cvxpy.log_sum_exp(exponents) + constant_term
                )

        # MSet variable constraints.
        for v in self._msets:
            xs, rhs = [], self._msets[v]
            for i, e in enumerate(rhs):
                matrix, coeffs, constant_term = e.specification(n)
                exponents = matrix @ variables + coeffs
                # xs.append(1/(i+1) * cvxpy.exp(cvxpy.sum(exponents)))
                # cvxpy.sum is not supported in Python2
                xs.append(1 / (i + 1) * cvxpy.exp(sum(exponents)))

            # constraints.append(variables[v.idx] >= cvxpy.sum(xs))
            # cvxpy.sum is not supported in Python2
            constraints.append(variables[v.idx] >= sum(xs))

        # Set variable constraints.
        for v in self._sets:
            xs, rhs = [], self._sets[v]
            for i, e in enumerate(rhs):
                matrix, coeffs, constant_term = e.specification(n)
                exponents = matrix @ variables + coeffs
                xs.append(sum(exponents))

            constraints.append(variables[v.idx] >= cvxpy.exp(sum(xs)))

        return constraints

    def _unfold_variables(self):
        """ Unfolds some more involved constructors.
        Note: this method might generate new variables in the system."""

        for v in self._seq_variables.copy():
            # note: since unfolding Seq might produce new sequence
            # variables, we cannot iterate over a set of dynamic
            # size. Instead, we iterate over its copy.
            v.register(self)

        for v in self._mset_variables:
            v.register(self)

        for v in self._ucyc_variables.copy():
            # same argument as for seq.
            v.register(self)

        for v in self._set_variables:
            v.register(self)

        for v in self._cyc_variables:
            v.register(self)

    def _run_solver(self, var, problem, params):
        """ Invokes the CVXPY solver."""

        if params.sys_type == Type.RATIONAL:
            solution = problem.solve(
                solver=params.solver,
                verbose=params.verbose,
                eps=params.eps,
                max_iters=params.max_iters,
            )
        else:
            solution = problem.solve(
                solver=params.solver,
                verbose=params.verbose,
                feastol=params.feastol,
                max_iters=params.max_iters,
                abstol=params.abstol,
                reltol=params.reltol,
            )

        # decorate system variables
        if var.value is not None:
            for idx, expr in enumerate(var.value):
                self._all_variables[idx].value = sympy.exp(expr).evalf()

        return solution

    def _inverse_synonym_variables(self):
        """ For each synonym variable, e.g. T = Seq(2z),
        includes Seq(2z) -> T in the specification dependency graph. """
        for v in self._equations:
            if self._is_synonym_variable(v):
                w = self._equations[v]._expressions[0]
                self._add_dependency(v, w)
                self._add_dependency(w, v)

    def _unreachable_nodes(self, target):
        """ Given the target variable, constructs the associated dependency
        graph and returns a list of nodes which cannot be reached from the
        target one."""

        self._build_dependency_graph()

        all = set(self._graph.nodes)
        ps = set(nx.shortest_path_length(self._graph, source=target.idx))

        return all - ps

    def _check_finite_tuner(self, target):
        """ Given the target variable, constructs the associated dependency
        graph and checks if all nodes are reachable from the target one."""

        self._build_dependency_graph()
        self._inverse_synonym_variables()

        unreachable = self._unreachable_nodes(target)
        return len(unreachable) == 0

    def _is_strongly_connected(self):
        """ Checks if the dependency graph is strongly connected. """

        self._build_dependency_graph()
        self._inverse_synonym_variables()

        d = dict(nx.shortest_path_length(self._graph))
        for v1 in self._equations:
            for v2 in self._equations:
                if v2.idx not in d[v1.idx]:
                    return False

        return True

    def _check_singular_tuner(self):
        """ Constructs the associated dependency graph and checks if it's
        strongly connected."""

        return self._is_strongly_connected()

    def run_tuner(self, t, params=None, method=Method.STRICT):
        """ Given the type variable and a set of tuning parameters, composes a
        (tuning) optimisation problem corresponding to an approximate sampler
        meant for structures of the given type. Variables are tuned so to
        achieve (in expectation) the marked variable values.

        Consider the following example:

          >>> sp = Specification()
          >>> z, u, M = Variable(1000), Variable(200), Variable()
          >>> sp.add(M, z + u * z * M + z * M **2)
          >>> params = Params(Type.ALGEBRAIC)
          >>> sp.run_tuner(M, params)

        Here, the variables z and u are marked with *absolute* values 1000 and
        200, respectively. The input type represents the type of Motzkin trees,
        i.e. plane unary-binary trees. Variable `z` marks their size, whereas
        `u` marks the occurrences of unary nodes. The tuning goal is to obtain
        specific values of z, u, and M, such that the induced branching
        probabilities lead to a sampler which generates Motzkin trees of size
        1000 with around 200 unary nodes (both in expectation).

        Respective variables (including type variables) are decorated with a
        proper 'value'. The method returns the CVXPY solution (i.e. the optimal
        value for the problem, or a string indicating why the problem could not
        be solved).

        By default, before the tuning procedure begins, certain sanity checks
        are performed, testing that the input specification matches necessary
        theoretical premises. It is possible to forcefully disable this
        behaviour by passing 'Method.FORCE' as the 'method' parameter."""

        assert len(self._tuned_variables) > 0, "No variables with tuning parameters."

        # get some default parameters if none given.
        params = self._init_params(params)

        # register unfoldable variables
        self._unfold_variables()

        if method == Method.STRICT:
            # check theoretical conditions imposed on finite-size tuning.
            if not self._check_finite_tuner(t):
                unreachable = self._unreachable_nodes(t)
                raise ValueError(
                    "Not all variables are reachable from "
                    "the target one. Please check for possible specification "
                    "errors or consider reformulating the specification. "
                    "Unreachable variables: " + str(unreachable) + ".\n"
                    "Note: indices correspond to variable names, see "
                    "repr(specification)."
                )

        n = self.discharged_variables
        variables = cvxpy.Variable(n)

        # compose the constraints
        constraints = self._compose_constraints(variables, n)

        # compose the objective
        obj = np.zeros(n, dtype="double")
        obj[t.idx] = 1.0

        for v in self._tuned_variables:
            obj[v.idx] = -v.tuning_param

        objective = cvxpy.Minimize(obj @ variables.T)
        problem = cvxpy.Problem(objective, constraints)
        return self._run_solver(variables, problem, params)

    def run_singular_tuner(self, z, params=None, method=Method.STRICT):
        """ Given a (size) variable and a set of tuning parameters, composes an
        optimisation problem corresponding to an approximate sampler meant for
        structures of the given type. Variables are tuned so to achieve (in
        expectation) the marked variable frequencies.

        Consider the following example:

          >>> sp = Specification()
          >>> z, u, M = Variable(), Variable(0.4), Variable()
          >>> sp.add(M, z + u * z * M + z * M **2)
          >>>
          >>> params = Params(Type.ALGEBRAIC)
          >>> sp.run_singular_tuner(z, M, params)

        Here, the variable u is marked with a *frequency* 0.4. The type M
        represents the type of Motzkin trees, i.e. unary-binary plane trees.
        Variable z marks their size, whereas u marks the occurrences of unary
        nodes. The tuning goal is to obtain specific values of `z`, `u`, and
        `M`, such that the induced branching probabilities lead to a sampler
        which generates, in expectation, Motzkin trees of infinite (i.e.
        unbounded) size and around 40% of unary nodes.

        Respective variables (including type variables) are decorated with a
        proper 'value'. The method returns the CVXPY solution (i.e. the optimal
        value for the problem, or a string indicating why the problem could not
        be solved).

        By default, before the tuning procedure begins, certain sanity checks
        are performed, testing that the input specification matches necessary
        theoretical premises. It is possible to forcefully disable this
        behaviour by passing 'Method.FORCE' as the 'method' parameter."""

        assert (
            z.tuning_param is None
        ), "Size parameter cannot be tuned in singular tuning."

        # get some default parameters if none given.
        params = self._init_params(params)

        # register unfoldable variables
        self._unfold_variables()

        if method == Method.STRICT:
            # check theoretical conditions imposed on singular tuning.
            if not self._check_singular_tuner():
                raise ValueError(
                    "Given specification has not a strongly connected "
                    "dependency graph. Please check the specification for "
                    "possible errors, or consider using finite-size "
                    "tuning for large values of the size parameter."
                )

        n = self.discharged_variables
        variables = cvxpy.Variable(n)

        # compose the constraints
        constraints = self._compose_constraints(variables, n)

        if params.sys_type == Type.RATIONAL:
            # for rational systems the optimisation problem becomes unbounded,
            # hence we need to artificially bound the vector norm.
            constraints.append(cvxpy.norm(variables, 2) <= params.norm)

        # compose the objective
        obj = np.zeros(n, dtype="double")
        obj[z.idx] = 1.0

        for v in self._tuned_variables:
            obj[v.idx] = v.tuning_param

        objective = cvxpy.Maximize(obj @ variables.T)
        problem = cvxpy.Problem(objective, constraints)
        return self._run_solver(variables, problem, params)

    def _variable_expressions(self, variable):
        if variable in self._msets:
            return Polynomial.sum(self._msets[variable])
        elif variable in self._sets:
            return Polynomial.sum(self._sets[variable])
        else:
            return self._equations[variable]

    def discrete_distribution(self, variable):
        """ Computes the discrete probability distribution for the given variable."""
        weights = [
            expr.weight / variable.value
            for expr in self._variable_expressions(variable)
        ]
        fractions = [Fraction(str(weight)).limit_denominator() for weight in weights]

        # make sure that the distribution sums up to one.
        (num, denum) = sum(fractions).as_integer_ratio()
        fractions[-1] += Fraction(denum - num, denum)

        for x in fractions:
            if x < Fraction(0, 1):
                raise ValueError(
                    "The distribution associated with the given variable contains negative components."
                )

        return fractions

    def approx_discrete_distribution(self, variable, precision=32, kernel="hellinger"):
        """ Computes an approximate discrete probability distribution for the given variable
        using the specified amount of precision bits and the kernel."""
        return get_optimal_probabilities(
            2 ** precision, self.discrete_distribution(variable), KERNELS[kernel]
        )

    def ddg(self, variable, precision=32, kernel="hellinger"):
        """ Computes a compact discrete distribution generating tree for the
        given variable. The result DDG uses an approximate discrete distribution.
        Note: Trivial distributions where all the probability mass is accumulated
        in a single point result in an empty list. """
        distribution = self.approx_discrete_distribution(variable, precision, kernel)
        for x in distribution:
            if x == Fraction(1, 1):
                return []

        enc, n, k = construct_sample_ky_encoding(distribution)
        return enc
