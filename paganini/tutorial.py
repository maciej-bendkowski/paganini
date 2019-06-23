"""
Tutorial
========

.. tip::
    Interactive environments like ``jupyter notebook`` are extremely helpful in
    code testing and experimenting. `Check them out! <https://jupyter.org>`_

.. note::
    Throughout the tutorial, it is assumed that at the beginning of the session,
    all the contents of the package `Paganini` have been imported:

        >>> from paganini import *

    Alternatively, in order to avoid polluting the global namespace, a
    synonym import can be used. In this case, all the functions should be
    referenced as sub-items of this namespace

        >>> import paganini as pg
        >>> spec = pg.Specification()

Introduction
------------

Consider the following example. Suppose that we are interested in designing an sampler for plane trees of unbounded degree (i.e. with an arbitrary number of children), specified as

    ``T = Z SEQ(T)``

where ``SEQ(T)`` stands for a (possibly empty) sequence of trees and ``Z`` marks the size of a node.
In Paganini, we write the following snippet defining the same combinatorial class:

    >>> spec = Specification()
    >>> z, T = Variable(), Variable()
    >>> spec.add(T, z * Seq(T))

Now, if we we want to construct a corresponding sampler, say *analytic (or Boltzmann) sampler*, we have to find a specific value of ``Z`` and use it to compute branching probabilities governing the random choices of our sampler (essentially the number of children for each of the constructed nodes). What value of ``Z`` should be choose if we are interested in large, uniform, and unbounded in size trees? With Paganini, this task amounts to invoking

    >>> spec.run_singular_tuner(z)

... and that's it! Paganini determines the corresponding value of `z` for us.
Once *tuned*, variables are decorated with appropriate numerical values:

    >>> z.value
    0.25
    >>> T.value
    0.5

Paganini allows its users to focus on the design of specifications, taking care of the rest.

Target expectation tuning
-------------------------

With the help of Paganini, users can demand the computation of tuning variables with specific, finite target expectations. Suppose that we are interested in designing an analytic sampler for Motzkin trees (i.e. plane unary-binary trees) however we would also like to demand that the outcome trees consists of around 1000 nodes, among which around 200 are  unary. To achieve this goal, we construct the following specification:

    >>> from paganini import *
    >>> spec = Specification()
    >>> z, u, M = Variable(1000), Variable(200), Variable()
    >>> spec.add(M, z + u * z * M + z * M ** 2)
    >>> spec.run_tuner(M)

Here ``z`` and ``u`` are two *marking variables* standing for the tree size and the number of unary nodes, respectively. Once we run the tuner, all three variables are decorated with respective numerical values, which the user can then use to compute respective branching probabilities. A sampler designed with such values is *guaranteed* to output Motzkin trees for which the expected size and mean number of unary nodes obey the design specification.


Examples
--------

Paganini is under constant development, supporting a growing class of so-called admissible constructors. Below you can find a handful of examples supported by Paganini. For more specifications, please visit our `tests` folder.


.. function:: Polya trees

    The specification is ``T = Z * MSET(T)``

    >>> spec = Specification()
    >>> z, T = Variable(), Variable()
    >>> spec.add(T, z * MSet(T))
    >>> spec.run_singular_tuner(z)
    >>> z.value
    0.338322112871298
    >>> T.value
    1.0

.. function:: Cyclic compositions

    A non-recursive specification ``C = CYC(Z * SEQ(Z))``
    with a dominant singularity ``z = 0.5``

    >>> spec = Specification()
    >>> z, C = Variable(), Variable()
    >>> spec.add(C, Cyc(z * Seq(z)))
    >>> spec.run_singular_tuner(z)
    >>> z.value
    0.499999992983969

.. function:: Unlabelled functional graphs.

    The specification is decined as a system of equations

    ``F = MSet(K)``

    ``K = CYC(U)``

    ``U = Z * MSet(U)``

    >>> spec = Specification()
    >>> z, F = Variable(), Variable()
    >>> K, U = Variable(), Variable()

    >>> spec.add(F, MSet(K))
    >>> spec.add(K, Cyc(U))
    >>> spec.add(U, z * MSet(U))

    >>> spec.run_singular_tuner(z)
    >>> z.value
    0.3383218568992077
"""
