"""
Paganini
========

Paganini is a lightweight python library for tuning of
multiparametric combinatorial systems.

All the necessary documentation can be found on-line on
https://paganini.readthedocs.io/

Use
    >>> help(paganini.tutorial)

to see some examples of code usage.
"""

__version__ = "1.2.1"

from .expressions import *
from .specification import *

import paganini.tutorial
import paganini.tests
