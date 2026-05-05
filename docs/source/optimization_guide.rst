KGBN Optimizer
==============

KGBN optimizes PBN rule probabilities against experimental measurements.

Overview
--------

Optimization requires:

- a ``ProbabilisticBN`` object,
- an experiment CSV or parsed experiment list,
- optional ``nodes_to_optimize`` to restrict which node probabilities can vary.

Basic Usage
-----------

.. code-block:: python

   import KGBN

   pbn = KGBN.load_pbn_from_file("network.txt")

   optimizer = KGBN.ParameterOptimizer(
       pbn,
       "experiments.csv",
       nodes_to_optimize=["Cas3"],
       verbose=False,
       normalize=False,
   )

   result = optimizer.optimize(method="differential_evolution")

   if result.success:
       optimized_pbn = optimizer.get_optimized_pbn(result)

Experimental Data
-----------------

Experiment CSV files describe perturbations and measurements:

.. code-block:: csv

   Experiments,Stimuli,Stimuli_efficacy,Inhibitors,Inhibitors_efficacy,Measured_nodes,Measured_values
   1,TGFa,1,TNFa,1,"NFkB,ERK,C8,Akt","0.7,0.88,0,1"
   2,TNFa,1,TGFa,1,"NFkB,ERK,C8,Akt","0.3,0.12,1,0"

Important columns:

- ``Stimuli``: nodes fixed to 1.
- ``Stimuli_efficacy``: optional efficacy values from 0 to 1.
- ``Inhibitors``: nodes fixed to 0.
- ``Inhibitors_efficacy``: optional efficacy values from 0 to 1.
- ``Measured_nodes``: measured node names or a formula.
- ``Measured_values``: target values for measured nodes or the formula.

Formula-Based Measurements
--------------------------

Use a formula when the objective should fit an aggregate phenotype-like score
instead of individual node values.

.. code-block:: python

   optimizer = KGBN.ParameterOptimizer(
       pbn,
       "experiments.csv",
       nodes_to_optimize=["nodeX"],
       Measured_formula="nodeA + nodeB - nodeC",
       normalize=True,
   )

Formulas support ``+``, ``-``, ``*``, ``/``, parentheses. Formula variables must match network node names.

Configuration
-------------

KGBN currently supports differential evolution and particle swarm optimization.
Pass method-specific parameters through ``config``:

.. code-block:: python

   config = {
       "seed": 9,
       "success_threshold": 0.005,
       "max_try": 3,
       "de_params": {
           "maxiter": 500,
           "popsize": 15,
           "workers": -1,
           "early_stopping": True,
       },
       "pso_params": {
           "n_particles": 30,
           "iters": 100,
           "options": {"c1": 0.5, "c2": 0.3, "w": 0.9},
       },
       "steady_state": {
           "method": "monte_carlo",
           "monte_carlo_params": {
               "n_runs": 10,
               "n_steps": 1000,
           },
       },
   }

   optimizer = KGBN.ParameterOptimizer(pbn, "experiments.csv", config=config)
   result = optimizer.optimize(method="differential_evolution")

Results And Evaluation
----------------------

The optimizer returns a SciPy-style result object with fields such as
``success``, ``message``, ``x``, ``fun``, ``history``, ``nfev``, and ``nit``.

Plot optimization history:

.. code-block:: python

   optimizer.plot_optimization_history(
       result,
       save_path="optimization_history.png",
       show_stagnation=True,
       log_scale=True,
   )

Evaluate the optimized model:

.. code-block:: python

   evaluator = KGBN.evaluate_optimization_result(
       result,
       optimizer,
       output_dir="evaluation_results",
       plot_residuals=True,
       save=True,
       detailed=True,
   )

See Also
--------

- :doc:`experiment_data`
- :doc:`parameter_optimizer`
- :doc:`result_evaluation`

References
----------

- Trairatphisan, P., et al. (2014). "optPBN: An Optimisation Toolbox for Probabilistic Boolean Networks." PLOS ONE 9(7): e98001.
