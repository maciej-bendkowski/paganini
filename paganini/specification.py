import cvxpy

from math import gcd
from paganini.expressions import *

class Seq(Variable):
    """ Sequence variables."""

    def __init__(self, expression, constraint = None):
        super(Seq,self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE # make sure its a type variable
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

        else: # constraint.operator == Operator.GEQ
            # note: Seq(expr)_{>= k} = expr^k + expr^{k+1} + ...
            #                        = expr^k (1 + expr^2 + expr^3 + ...)
            #                        = expr^k Seq(expr).

            seq = Seq(self.inner_expressions)
            spec.add(seq, 1 + self.inner_expressions * seq)
            spec.add(self, self.inner_expressions ** self.constraint.value * seq)

class MSet(Variable):
    """ MSet variables."""

    def __init__(self, expression):
        super(MSet,self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE # make sure its a type variable

    def register(self, spec):
        """ Unfolds the MSet definition and registers it in the given system."""
        spec._diagonal_variable(self, 1)

class Cyc(Variable):
    """ Cyc variables."""

    def __init__(self, expression):
        super(Cyc,self).__init__()

        expression = Polynomial.cast(expression)
        self.inner_expressions = expression
        self.type = VariableType.TYPE # make sure its a type variable

    def register(self, spec):
        """ Unfolds the Cyc definition and registers it in the given system."""
        spec._diagonal_variable(self, 1)

class Operator(Enum):
    """ Enumeration of supported constraint signs."""
    LEQ       = 1 # less or equal
    GEQ       = 2 # greater or equal
    UNBOUNDED = 3 # unbounded operator

class Constraint:
    """ Supported constraints for classes such as SEQ."""
    def __init__(self, operator, value):
        self.operator = operator
        self.value    = value

    @staticmethod
    def normalise(constraint = None):

        if constraint is None:
            return Constraint(Operator.UNBOUNDED, 0)
        else:
            return constraint

def leq(n):
    """ Creates a less or equal constraint for the given input."""
    assert n >= 0, "Negative constraints are not supported."
    return Constraint(Operator.LEQ, n)

def geq(n):
    """ Creates a greater or equal constraint for the given input."""
    assert n >= 0, "Negative constraints are not supported."
    return Constraint(Operator.GEQ, n)

class Type(Enum):
    """ Enumeration of supported system types."""
    ALGEBRAIC = 1
    RATIONAL  = 2

class Params:
    """ CVXPY solver parameters initalised with some defaults."""

    def __init__(self, sys_type):
        self.verbose   = True
        if sys_type == Type.RATIONAL:
            self.sys_type  = Type.RATIONAL
            self.solver    = cvxpy.SCS
            self.max_iters = 2500
            self.eps       = 1.e-20
            self.norm      = 40
        else:
            self.sys_type  = Type.ALGEBRAIC
            self.solver    = cvxpy.ECOS
            self.max_iters = 100
            self.feastol   = 1.e-20

def phi(n):
    """ Euler's totient function."""
    assert n >= 0, 'Negative integer.'

    out = 0
    for i in range(1, n + 1):
        if gcd(n, i) == 1:
            out += 1

    return out

class Specification:
    """ Symbolic system specifications."""

    def __init__(self, truncate = 10):
        """ Creates a new specification. The optional `truncate` parameter
        controls the truncation threshold for infinite series, intrinsic to
        some more involved constructions, such as multisets or cycles."""

        self._index_counter = 0
        self._equations = {}

        self._all_variables   = {}
        self._tuned_variables = set()

        self._seq_variables   = set()
        self._mset_variables  = set()
        self._cyc_variables   = set()

        # diagonals.
        self._diag     = {}
        self._msets    = {}
        self._cycs     = {}

        self._truncate = truncate

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
            return # nothing to do.

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

        elif isinstance(v, Cyc):
            self._cyc_variables.add(v)
            self._register_expressions(v.inner_expressions)

    def _register_expression(self, expression):
        assert isinstance(expression, Expr), 'Expected expression.'
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

    def check_type(self):
        """ Checks if the system is algebraic or rational.
        Note: the current method is heuristic."""

        if len(self._seq_variables) > 0 or\
                len(self._mset_variables) > 0 or\
                len(self._cyc_variables) > 0:
            return Type.ALGEBRAIC

        for expressions in self._equations.values():
            for exp in expressions:
                if isinstance(exp, Expr):
                    for v, e in exp.variables.items():
                        if v.is_type_variable and e > 1:
                            return Type.ALGEBRAIC

        return Type.RATIONAL

    def _init_params(self, params = None):
        if params is None:
            # some defaults
            sys_type = self.check_type()
            return Params(sys_type)
        else:
            return params

    def _diagonal_expr(self, expr, d = 1):
        """ Extends self._diag_variable to expressions."""
        assert isinstance(expr, Expr), 'Expression required.'

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

    def _diagonal_variable(self, var, d = 1):
        """ Given T(z) creates its kth diagonal, i.e. T(z^k)."""

        assert d > 0, 'Non-positive diagonal parameter.'
        assert var.is_type_variable, 'Requested diagonal of non-type variable.'

        if var not in self._diag:
            self._diag[var] = {} # degree -> expression

        # check the variable cache.
        if d in self._diag[var]:
            return self._diag[var][d]

        # trivial diagonal case.
        if d == 1 and var in self._equations:
            self._diag[var][d] = var
            return var

        var_d = self._type_variable() if d > 1 else var
        self._diag[var][d] = var_d # memorise var[d]

        if var in self._equations:
            monomials = self._equations[var]
            self.add(var_d, [self._diagonal_expr(e, d) for e in monomials])
            return var_d

        elif var in self._mset_variables:
            self._msets[var_d] = []
            for k in range(d, self._truncate + 1, d):
                self._msets[var_d].append(Polynomial([self._diagonal_expr(e, k)\
                        for e in var.inner_expressions]))

            return var_d

        elif var in self._cyc_variables:
            self._cycs[var_d] = []
            for k in range(d, self._truncate + 1, d):

                # create a sequence version of the expressions
                seq_k = Seq([self._diagonal_expr(e, k) for e in
                            var.inner_expressions])

                self._register_variable(seq_k)
                seq_k.register(self)

                self._cycs[var_d].append(Polynomial([seq_k]))

            return var_d

    def add(self, var, expression):
        """ Includes the given definition in the specification."""
        expression = Polynomial.cast(expression)
        var.type = VariableType.TYPE # make var a type variable
        self._equations[var] = expression

        # register variables in the system.
        self._register_variable(var)
        self._register_expressions(expression)

    def _compose_constraints(self, variables, n):
        """ Composes optimisation constraints."""

        constraints = []
        for v in self._equations:
            rhs = self._equations[v] # right-hand side.
            matrix, coeffs, constant_term = rhs.specification(n)

            exponents = matrix * variables + coeffs
            constraints.append(variables[v.idx] >=
                cvxpy.log_sum_exp(exponents) + constant_term)

        # MSet variable constraints.
        for v in self._msets:
            xs, rhs = [],  self._msets[v]
            for i, e in enumerate(rhs):
                matrix, coeffs, constant_term = e.specification(n)
                exponents = matrix * variables + coeffs
                xs.append(1/(i+1) * cvxpy.exp(cvxpy.sum(exponents)))

            constraints.append(variables[v.idx] >= cvxpy.sum(xs))

        # Cyc variable constraints.
        for v in self._cycs:
            xs, rhs = [],  self._cycs[v]
            for i, e in enumerate(rhs):
                matrix, coeffs, constant_term = e.specification(n)
                exponents = matrix * variables + coeffs
                xs.append(phi(i+1)/(i+1) * exponents)

            constraints.append(variables[v.idx] >= cvxpy.sum(xs))

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

        for v in self._cyc_variables:
            v.register(self)

    def _run_solver(self, var, problem, params):
        """ Invokes the CVXPY solver."""

        if params.sys_type == Type.RATIONAL:
            solution = problem.solve(solver = params.solver, verbose =
                    params.verbose, eps = params.eps, max_iters =
                    params.max_iters)
        else:
            solution = problem.solve(solver = params.solver, verbose =
                    params.verbose, feastol = params.feastol, max_iters =
                    params.max_iters)

        # decorate system variables
        for idx, expr in enumerate(var.value):
            self._all_variables[idx].value = sympy.exp(expr).evalf()

        return solution

    def run_tuner(self, t, params = None):
        """ Given the type variable and a set of tuning parameters, composes a
        (tuning) optimisation problem corresponding to an approximate sampler
        meant for structures of the given type. Variables are tuned so to
        achieve (in expectation) the marked variable values.

        Consider the following example:

        sp = Specification()
        z, u, M = Variable(1000), Variable(200), Variable()
        sp.add(M, z + u * z * M + z * M **2)
        params = Params(Type.ALGEBRAIC)
        sp.run_tuner(M, params)

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
        be solved)."""

        assert len(self._tuned_variables) > 0,\
                'No variables with tuning parameters.'

        # get some default parameters if none given.
        params = self._init_params(params)

        # register unfoldable variables
        self._unfold_variables()

        n = self.discharged_variables
        variables = cvxpy.Variable(n)

        # compose the constraints
        constraints = self._compose_constraints(variables, n)

        # compose the objective
        obj = np.zeros(n)
        obj[t.idx] = 1.0

        for v in self._tuned_variables:
            obj[v.idx] = - v.tuning_param

        objective = cvxpy.Minimize(obj * variables)
        problem   = cvxpy.Problem(objective, constraints)
        return self._run_solver(variables, problem, params)

    def run_singular_tuner(self, z, params = None):
        """ Given a (size) variable and a set of tuning parameters, composes an
        optimisation problem corresponding to an approximate sampler meant for
        structures of the given type. Variables are tuned so to achieve (in
        expectation) the marked variable frequencies.

        Consider the following example:

        sp = Specification()
        z, u, M = Variable(), Variable(0.4), Variable()
        sp.add(M, z + u * z * M + z * M **2)

        params = Params(Type.ALGEBRAIC)
        sp.run_singular_tuner(z, params)

        Here, the variable u is marked with a *frequency* 0.4.  The type M
        represents the type of Motzkin trees, i.e. unary-binary plane trees.
        Variable z marks their size, whereas u marks the occurrences of unary
        nodes. The tuning goal is to obtain specific values of z, u, and M, such
        that the induced branching probabilities lead to a sampler which
        generates, in expectation, Motzkin trees of infinite (i.e. unbounded)
        size and around 40% of unary nodes.

        Respective variables (including type variables) are decorated with a
        proper 'value'. The method returns the CVXPY solution (i.e. the optimal
        value for the problem, or a string indicating why the problem could not
        be solved)."""

        # get some default parameters if none given.
        params = self._init_params(params)

        # register unfoldable variables
        self._unfold_variables()

        n = self.discharged_variables
        variables = cvxpy.Variable(n)

        # compose the constraints
        constraints = self._compose_constraints(variables, n)

        if params.sys_type == Type.RATIONAL:
            # for rational systems the optimisation problem becomes unbounded,
            # hence we need to artificially bound the vector norm.
            constraints.append(cvxpy.norm(variables, 2) <= params.norm)

        # compose the objective
        obj = np.zeros(n)
        obj[z.idx] = 1.0

        for v in self._tuned_variables:
            obj[v.idx] = v.tuning_param

        objective = cvxpy.Maximize(obj * variables)
        problem   = cvxpy.Problem(objective, constraints)
        return self._run_solver(variables, problem, params)
