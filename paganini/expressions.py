import sympy
import numpy as np

from enum import Enum
from scipy import sparse
from collections import Counter

__all__ = ('Expr', 'VariableType', 'Variable', 'Polynomial')

class Expr(object):
    """ Algebraic expressions (multivariate monomials) in form of
        :math:`c x_1^{k_1} x_2^{k_2} \cdots x_m^{k_m}.`"""

    def __init__(self, coeff = 1, variables = Counter()):
        self.variables = variables
        self.coeff = coeff

    @staticmethod
    def cast(other):
        """ Casts its input to an expression."""

        if not isinstance(other, Expr):
            return Expr(other)
        else:
            return other

    def __add__(self, other):
        """ Expression addition."""
        return Polynomial(self) + other # lift to polynomials.

    __radd__ = __add__ # make addition commute again

    def __sub__(self, other):
        """ Expression subtraction."""
        return Polynomial(self) - other # lift to polynomials.

    def __mul__(self, other):
        """ Multiplication of algebraic expressions."""
        if isinstance(other, Polynomial):
            return Polynomial(self) * other

        other = Expr.cast(other)

        # a^b * a^c = a^(b + c)
        x, y = Counter(self.variables), Counter(other.variables)
        return Expr(self.coeff * other.coeff, x + y)

    __rmul__ = __mul__ # make multiplication commute again.

    def __pow__(self, n):
        """ Expression exponentiation."""
        assert n > 0, 'Non-positive exponent.'

        xs = Counter(self.variables)
        for v in self.variables:
            xs[v] *= n # (a^b)^c = a^(b * c)

        return Expr(self.coeff ** n, xs)

    def __repr__(self):
        monic_monomial = ' '.join([
                variable.__repr__() + "^" + str(self.variables[variable])
                for variable in self.variables])

        if self.coeff == 1 and len(monic_monomial) > 0:
            return monic_monomial
        else:
            return str(self.coeff) + " " + monic_monomial

    def related(self, other):
        """ Checks if the two expressions are related, i.e. have
        the same exponents and perpahs a different coefficient."""
        other = Expr.cast(other)
        return self.variables == other.variables

    @property
    def is_constant(self):
        """ True iff the expression represents a constant."""
        return len(self.variables) == 0

class VariableType(Enum):
    PLAIN = 1 # regular, plain variables, e.g. Z.
    TYPE  = 2 # variables corresponding to some types, i.e. having definitions.

class Variable(Expr):
    """ Symbolic variables."""

    def __init__(self, tuning_param = None):
        super(Variable, self).__init__(1, Counter())

        self.variables[self] = 1
        self.type  = VariableType.PLAIN
        self.tuning_param = tuning_param

        self.idx   = None
        self.value = None

    @property
    def is_type_variable(self):
        """ True iff the variable represents a type variable.
        In other words, if it admits a defining equation."""
        return self.type == VariableType.TYPE

    def __repr__(self):
        return "var" + str(self.idx)

    __str__ = __repr__

    def set_expectation(self, tuning_param):
        self.tuning_param = tuning_param

class Polynomial:
    """ Polynomials of multivariate algebraic expressions."""

    def __init__(self, expressions):
        if isinstance(expressions, Expr):
            expressions = [expressions]

        self._expressions = expressions

    @staticmethod
    def cast(other):
        """ Casts its input to a polynomial."""
        if isinstance(other, (int,float)):
            return Polynomial([Expr(other)])

        elif not isinstance(other, Polynomial):
            return Polynomial(other)

        else:
            return other

    @staticmethod
    def simplify(polynomial):
        """ Simplifies the given polynomial."""

        equiv_classes = []
        n = len(polynomial)
        visited = [False] * n

        # group expressions
        for i in range(0, n):
            if not visited[i]:
                visited[i] = True
                equiv_class = [polynomial[i]]

                for j in range(i + 1, n):
                    if polynomial[i].related(polynomial[j]):
                        equiv_class.append(polynomial[j])
                        visited[j] = True

                equiv_classes.append(equiv_class)

        # collect coefficients
        simpl_expressions = []
        for eqv in equiv_classes:

            coeff = 0
            for expr in eqv:
                coeff += expr.coeff

            if coeff != 0: # ignore vacuous terms
                expr = Expr(coeff, eqv[0].variables)
                simpl_expressions.append(expr)

        return Polynomial(simpl_expressions)

    @staticmethod
    def sum(series):
        """ Evaluates the sum of the given series."""

        if len(series) > 0:
            p = series[0]
            for i in range(1, len(series)):
                p += series[i]

            return p
        else:
            return Polynomial(Expr(0))

    def __add__(self, other):
        """ Polynomial addition."""
        other = Polynomial.cast(other)
        return Polynomial.simplify(self._expressions + other._expressions)

    __radd__ = __add__ # make addition commute again

    def __sub__(self, other):
        """ Polynomial subtraction."""
        if isinstance(other, (int,float)):
            return self + (-other)

        other = Polynomial.cast(other)
        xs = [-1 * e for e in other._expressions]
        return self + Polynomial(xs)

    def __mul__(self, other):
        """ Naive polynomial multiplication."""
        other = Polynomial.cast(other)

        outcome = [] # naive but works
        for a in self._expressions:
            for b in other._expressions:
                outcome.append(a * b)

        return Polynomial.simplify(outcome)

    __rmul__ = __mul__ # make multiplication commute again

    def __pow__(self, n):
        """ Naive polynomial exponentiation."""
        assert n >= 0, 'Non-positive exponent.'

        if n == 0:
            return Polynomial(Expr(1))

        if n == 1:
            return Polynomial(self._expressions)

        if n % 2 == 1:
            return self * self ** (n - 1)
        else:
            other = self ** (n >> 1)
            return other * other

    def __iter__(self):
        return iter(self._expressions)

    def __repr__(self):
        return ' + '.join([
            expression.__repr__()
            for expression in self._expressions
        ])

    def is_one(self):
        """ Checks if the polynomial represents a constant one."""
        return len(self._expressions) == 1\
                and self._expressions[0].is_constant\
                and self._expressions[0].coeff == 1

    def specification(self, no_variables):
        """ Composes a sparse matrix specification of the polynomial. Requires
        as input a number dictating the number of columns of the constructed
        matrix (usually the number of variables in the corresponding
        optimisation problem).

        Its output is a tuple consisting of:

        (1) a sparse matrix representing the polynomial,
        (2) a vector of logarithms of monomial coefficients,
        (3) a (collective) constant term representing constant monomials.

        The matrix represents expoenents of respective variables."""

        rows = 0 # row counter
        row, col, data =  [], [], []
        constant_expr = 0
        coeffs = []

        for exp in self:
            if isinstance(exp, Expr):
                if exp.coeff <= 0 and not exp.is_constant:
                    raise ValueError('Non-positive monomial coefficient.')

                if exp.coeff > 0:
                    #coeffs.append(sympy.log(exp.coeff))
                    coeffs.append(np.float64((sympy.log(exp.coeff)).evalf()))
                else:
                    constant_expr += exp.coeff

                for (v, e) in exp.variables.items():
                    row.append(rows)
                    col.append(v.idx)
                    data.append(e)
                rows += 1
            else:
                constant_expr += exp # constant

        # create a sparse representation of the polynomial,
        # together with logarithms of respective monomial coefficients and
        # the collected constant term (unaltered).
        return (sparse.csr_matrix((np.array(data),
            (np.array(row),np.array(col))), shape=(rows, no_variables),
            dtype='double'),
            np.array(coeffs), constant_expr)
