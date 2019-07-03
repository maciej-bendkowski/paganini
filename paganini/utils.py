from sympy import gcd

__all__ = ('phi')

def phi(n):
    """ Euler's totient function."""
    assert n >= 0, 'Negative integer.'

    out = 0
    for i in range(1, n + 1):
        if gcd(n, i) == 1:
            out += 1

    return out
