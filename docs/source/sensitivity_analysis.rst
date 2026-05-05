Sensitivity Analysis
====================

The ``KGBN.sensitivity_analysis`` module provides three helpers for analyzing
PBNs:

- ``sensitivity_analysis``: Morris or Sobol sensitivity analysis against
  experimental measurements.
- ``influence_analysis``: node-to-node influence analysis from steady-state
  samples.
- ``mse_sensitivity_analysis``: MSE change after perturbing probabilities in an
  optimized PBN.

Install optional dependencies with ``pip install "KGBN[sensitivity]"`` or
``pip install "KGBN[all]"``.

.. automodule:: KGBN.sensitivity_analysis
   :members:
   :undoc-members:
   :show-inheritance:

Inputs
------

Sensitivity functions expect a ``ProbabilisticBN`` object. Convert a BN to a
PBN before analysis:

.. code-block:: python

   import KGBN

   bn = KGBN.load_network("model.txt", network_type="bn")
   pbn_string, nodes_to_optimize = KGBN.BN2PBN(bn, prob=0.8)
   pbn = KGBN.load_network(pbn_string, network_type="pbn")

For ``sensitivity_analysis`` and ``mse_sensitivity_analysis``, pass either an
experiment CSV path or a list of experiment dictionaries. See
:doc:`experiment_data` for the CSV format.

Functions
---------

sensitivity_analysis
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: KGBN.sensitivity_analysis.sensitivity_analysis

.. code-block:: python

   from KGBN.sensitivity_analysis import sensitivity_analysis

   results = sensitivity_analysis(
       pbn,
       "experiments.csv",
       config={
           "sensitivity_method": "morris",  # or "sobol"
           "morris_trajectories": 20,
           "sobol_samples": 512,
           "parallel": True,
           "n_workers": -1,
           "steady_state": {
               "method": "tsmc",
               "tsmc_params": {
                   "epsilon": 0.01,
                   "max_iterations": 1000,
               },
           },
       },
       top_n=10,
   )

   print(results["sensitivity_df"])
   print(results["top_nodes"])

influence_analysis
~~~~~~~~~~~~~~~~~~

.. autofunction:: KGBN.sensitivity_analysis.influence_analysis

.. code-block:: python

   from KGBN.sensitivity_analysis import influence_analysis, plot_influence_results

   results = influence_analysis(
       pbn,
       config={
           "method": "tsmc",
           "samples": 10000,
           "burn_in": 2000,
           "thin": 5,
           "epsilon": 0.001,
       },
       top_n=10,
   )

   plot_influence_results(
       results["influence_df"],
       results["sensitivity_df"],
       top_n=20,
       save_path="influence_analysis.png",
   )

mse_sensitivity_analysis
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: KGBN.sensitivity_analysis.mse_sensitivity_analysis

.. code-block:: python

   from KGBN.sensitivity_analysis import (
       mse_sensitivity_analysis,
       plot_mse_sensitivity_results,
   )

   results = mse_sensitivity_analysis(
       optimized_pbn,
       "experiments.csv",
       config={
           "perturbation_magnitude": 0.3,
           "steady_state": {
               "method": "tsmc",
               "tsmc_params": {
                   "epsilon": 0.01,
                   "max_iterations": 1000,
               },
           },
       },
       top_n=10,
   )

   plot_mse_sensitivity_results(
       results["sensitivity_df"],
       results["baseline_mse"],
       top_n=20,
       save_path="mse_sensitivity.png",
   )

Formula Measurements
--------------------

If the objective uses a formula instead of direct measured nodes, pass
``Measured_formula`` and optionally enable min-max normalization:

.. code-block:: python

   results = sensitivity_analysis(
       pbn,
       "experiments.csv",
       Measured_formula="N1 + N2 - N3",
       normalize=True,
   )

See Also
--------

- :doc:`parameter_optimizer`
- :doc:`experiment_data`
- :doc:`steady_state_guide`
