Installation
============

`Paganini` is available as a `Python` package for both `Python2` and `Python3`.

.. warning::
    We strongly recommend using `Python3`, and
    we do not guarantee the required decimal precision for `Python2`.
    Moreover, the support for Python2 drops in 2020, and as a consequence, some
    of the packages that we use are not anymore maintained, and the usage is
    at your own risk.
    In particular, two of our tests fail on Python2 for the reasons of numerical
    precision.

.. note::
    We assume that the user is familiar with `Python` and its package manager
    `pip`. If you are new to `Python`, please visit the official installation
    webpage `<https://www.python.org/downloads/>`_. For new versions of `Python`,
    `pip` is already pre-installed. If you don't have it, check
    `<https://pip.pypa.io/en/stable/installing/>`_


The latest release of `Paganini` can be installed using `pip`:

    >>> pip install paganini

If you want to update to the recent version of paganini, use

    >>> pip install --upgrade paganini

Installation from sources
-------------------------

In order to install from sources, you need `git`, or you can download and
install the code manually from `github`.

::

    git clone git://github.com/maciej-bendkowski/paganini.git
    cd paganini
    python3 setup.py install

Testing
-------

In order to verify that paganini works, follow the steps of the tutorial

    >>> import paganini
    >>> help(paganini.tutorial)

