# Paganini
Paganini is a lightweight python library meant to help with the design of
combinatorial samplers. Given a combinatorial specification, expressed using a
domain-specific language closely resembling Flajolet and Sedgewick's *symbolic
method*, Paganini gives its users some additional control over the distribution
of structures constructed using the designed samplers.

### Citing Paganini
If you use Paganini or its components for published work, we encourage you to
cite the accompanying papers:

*Maciej Bendkowski, Olivier Bodini, Sergey Dovgal*

- [Tuning as convex optimisation: a polynomial tuner for multi-parametric combinatorial samplers](https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/tuning-as-convex-optimisation-a-polynomial-tuner-for-multiparametric-combinatorial-samplers/2A79BF69E8CC80011E6AADF81BFC0291)

- [Polynomial tuning of multiparametric combinatorial samplers](https://epubs.siam.org/doi/10.1137/1.9781611975062.9)

### Introduction

Consider the following example. Suppose that we are interested in designing a
sampler for plane trees of unbounded degree (i.e. with an arbitrary number of
children), specified as

```
T = Z SEQ(T)
```
where `SEQ(T)` stands for a (possibly empty) sequence of trees and `Z` marks the size of a node.
In Paganini, we write the following snippet defining the same combinatorial class:

```python
sp = Specification()
z, T = Variable(), Variable()
sp.add(T, z * Seq(T))
```

Now, if we we want to construct a corresponding sampler, say *analytic (or
Boltzmann) sampler*, we have to find a specific value of `Z` and use it to
compute branching probabilities governing the random choices of our sampler
(essentially the number of children for each of the constructed nodes). What
value of `Z` should be choose if we are interested in large, uniform, and
unbounded in size trees? With Paganini, this task amounts to invoking

```python
sp.run_singular_tuner(z)
```

... and that's it! Paganini determines the corresponding value of `z` for us.
Once *tuned*, variables are decorated with appropriate numerical values:

```python
z.value // 0.25
T.value // 0.5
```
Paganini allows its users to focus on the design of specifications, taking care of the rest.

### Target expectation tuning
With the help of Paganini, users can demand the computation of tuning variables
with specific, finite target expectations. Suppose that we are interested in
designing an analytic sampler for Motzkin trees (i.e. plane unary-binary trees)
however we would also like to demand that the outcome trees consists of around
1000 nodes, among which around 200 are unary. To achieve this goal, we construct
the following specification:

```python
sp = Specification()
z, u, M = Variable(1000), Variable(200), Variable()
sp.add(M, z + u * z * M + z * M ** 2)
sp.run_tuner(M)
```
Here `z` and `u` are two *marking variables* standing for the tree size
and the number of unary nodes, respectively. Once we run the tuner, all three variables
are decorated with respective numerical values, which the user can then use to compute
respective branching probabilities. A sampler designed with such values is *guaranteed*
to output Motzkin trees for which the expected size and mean number of unary nodes obey
the design specification.

### Examples
Paganini is under constant development, supporting a growing class of so-called
admissible constructors. For more specifications, please visit our

[Online tutorial](https://paganini.readthedocs.io/en/latest/tutorial.html)

### Installation
Paganini is available as a pip package.
```
pip install paganini
```

### References
Paganini relies on published work of numerous excellent authors. Below, you can
find a short (and definitely inexhaustive) list of papers on the subject:

- [P. Flajolet, R. Sedgewick: Analytic Combinatorics](http://algo.inria.fr/flajolet/Publications/book.pdf)
- [P. Duchon, P. Flajolet, G. Louchard. G. Schaeffer: Boltzmann Samplers for
   the random generation of combinatorial structures](http://algo.inria.fr/flajolet/Publications/DuFlLoSc04.pdf)
- [C. Pivoteau, B. Salvy, M. Soria: Algorithms for Combinatorial Systems:
   Well-Founded Systems and Newton Iterations](https://www.sciencedirect.com/science/article/pii/S0097316512000908)
- [O.Bodini, J. Lumbroso, N. Rolin: Analytic samplers and the combinatorial rejection method](https://dl.acm.org/citation.cfm?id=2790220&dl=ACM&coll=DL)

If you are interested in the practical design of analytic samplers, we encourage
you to check out the related [Boltzmann
Brain](https://github.com/maciej-bendkowski/boltzmann-brain) software.
