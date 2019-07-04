from sympy import gcd

__all__ = ('phi', 'partition_sequences')

def phi(n):
    """ Euler's totient function."""
    assert n >= 0, 'Negative integer.'

    out = 0
    for i in range(1, n + 1):
        if gcd(n, i) == 1:
            out += 1

    return out

def partition_sequences(k):
    """ Generates a set of partition_sequences of size :math:`k`, i.e.
    sequences :math:`(n_i)` such that :math:`\sum_{i=1}^k i n_i = k`."""
    assert k >= 0, 'Negative integer.'

    def f(k, c, i):

        if k == 0:
            yield [0] * c

        elif i == k:
            yield [1] + [0] * (c - 1)

        elif i < k:
            for n in range(0, k // i + 1):
                for ps in f(k - (i * n), c - 1, i + 1):
                    yield [n] + ps

    return f(k, k, 1)
