"""
Paganini
========

Paganini is a lightweight python library for tuning of
multiparametric combinatorial systems.

How to use the documentation
----------------------------

The docstring assumes that `paganini` is imported as `pg`::

  >>> import paganini as pg

General-purpose help is available in the `docs` subpackage::

  >>> help(pg.docs)

All the necessary documentation can be found on-line on
https://paganini.readthedocs.io/
"""

__version__ = "1.1.1"

from .expressions import *
from .specification import *
