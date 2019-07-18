import os
from setuptools import setup, find_packages

setup(
    name = "paganini",
    version = "1.2.1",

    author = "Maciej Bendkowski, Sergey Dovgal",
    author_email = "maciej.bendkowski@tcs.uj.edu.pl, vit.north@gmail.com",
    description = "Multiparametric tuner for combinatorial specifications",

    license = "BSD3",
    url = "https://github.com/maciej-bendkowski/paganini",
    install_requires = ['numpy', 'sympy', 'cvxpy', 'scipy'],
    packages = find_packages(),

    # unit tests
    test_suite = 'paganini.tests'
)
