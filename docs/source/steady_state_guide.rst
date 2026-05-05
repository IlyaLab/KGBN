Steady State Calculation Guide
==============================

This guide summarizes how to compute steady-state behavior for Boolean
Networks (BNs) and Probabilistic Boolean Networks (PBNs) in KGBN.

Overview
--------

Use ``KGBN.SteadyStateCalculator`` after loading a network. The calculator
selects deterministic attractor finding for BNs and supports Monte Carlo or
two-state Markov chain (TSMC) estimation for PBNs.

Quick Start
-----------

.. code-block:: python

   import KGBN

   pbn_string = """
   N1 = N1, 1
   N2 = N2, 1
   N3 = N1, 0.6
   N3 = N1 & !N2, 0.4
   """

   pbn = KGBN.load_network(pbn_string, initial_state=[1, 1, 0])
   calc = KGBN.SteadyStateCalculator(pbn)

   steady_state = calc.compute_steady_state(
       method="monte_carlo",
       n_runs=20,
       n_steps=10000,
       p_noise=0.05,
       seed=9,
   )

   print(steady_state)

Methods
-------

Monte Carlo
~~~~~~~~~~~

Monte Carlo runs independent simulations and averages the second half of each
trajectory.

Common parameters:

- ``n_runs``: number of independent runs.
- ``n_steps``: simulation steps per run.
- ``p_noise``: probability of random state flips.
- ``seed``: random seed.

.. code-block:: python

   steady_state = calc.compute_steady_state(
       method="monte_carlo",
       n_runs=20,
       n_steps=10000,
       p_noise=0.05,
       seed=9,
   )

TSMC
~~~~

TSMC estimates steady state with transition-probability criteria and can be
faster for large PBNs.

Common parameters:

- ``epsilon``: convergence threshold.
- ``r``: accuracy range.
- ``s``: probability of accuracy.
- ``p_mir``: perturbation probability.
- ``initial_nsteps``: initial simulation length.
- ``max_iterations``: maximum convergence iterations.
- ``freeze_constant``: whether self-loop nodes should stay fixed.

.. code-block:: python

   steady_state = calc.compute_steady_state(
       method="tsmc",
       epsilon=0.001,
       r=0.01,
       p_mir=0.01,
       initial_nsteps=100,
       max_iterations=5000,
       freeze_constant=False,
       seed=9,
   )

Deterministic Attractors
~~~~~~~~~~~~~~~~~~~~~~~~

For deterministic BNs, KGBN finds fixed points and cyclic attractors.

.. code-block:: python

   bn = KGBN.load_network("network.txt", network_type="bn")
   calc = KGBN.SteadyStateCalculator(bn)

   attractors = calc.compute_steady_state(
       method="deterministic",
       n_runs=100,
       n_steps=1000,
       verbose=True,
       seed=9,
   )

   print(attractors["fixed_points"])
   print(attractors["cyclic_attractors"])

Experimental Conditions
-----------------------

Use experimental conditions to fix stimuli to 1 and inhibitors to 0 before
steady-state calculation. Partial efficacy values are supported for Monte Carlo
simulation.

.. code-block:: python

   calc.set_experimental_conditions(
       stimuli=["N1"],
       stimuli_efficacy=[0.8],
       inhibitors=["N2"],
       inhibitors_efficacy=[0.5],
   )

   perturbed_state = calc.compute_steady_state(
       method="monte_carlo",
       n_runs=20,
       n_steps=10000,
       p_noise=0.05,
       seed=9,
   )

   calc.reset_network_conditions()

See Also
--------

- :doc:`basic_simulation`
- :doc:`optimization_guide`
- :doc:`steady_state`
