# KGBN: Knowledge Graph-Augmented Boolean Network Modeling

KGBN is a Python toolkit for augmenting and optimizing logical gene regulatory network models with knowledge-graph-derived regulatory hypotheses. It supports Boolean networks (BNs), probabilistic Boolean networks (PBNs), KG-based model extension, phenotype scoring, parameter optimization, simulation, evaluation, and visualization.

The package follows the workflow described in the KGBN manuscript:

1. **Simulate Boolean networks and probabilistic Boolean networks**
2. **Extend the model via knowledge graphs**
3. **Represent alternative regulatory hypotheses as PBN rules**
4. **Optimize PBN rule probabilities using experimental data**
5. **Simulate and evaluate the optimized model**

Documentation is available at [https://ilyalab.github.io/KGBN/](https://ilyalab.github.io/KGBN/).

## Installation

For development:

```bash
pip install -e .
```

## Quick Start

```python
import KGBN

network_string = """
MYC = MYC
CDKN2A = !MYC | TP53
TP53 = !CDKN2A
"""

bn = KGBN.load_network(network_string, initial_state={"MYC": 1, "CDKN2A": 0, "TP53": 0})
trajectory = bn.update(iterations=5)

calc = KGBN.SteadyStateCalculator(bn)
calc.set_experimental_conditions(stimuli=['MYC'])
steady_state = calc.compute_steady_state(method="monte_carlo", n_runs=3, n_steps=100)
```

## Knowledge Graph Extension

```python
genes = ["TP53", "MYC", "CDKN2A", "MAPK1"]
kg_rules, relations = KGBN.load_signor_network(
    genes,
    joiner="inhibitor_wins",
    score_cutoff=0.5,
)
kg_bn = KGBN.load_network(kg_rules)
print(kg_rules)
```

## PBN Optimization

```python
pbn_string = KGBN.merge_networks([bn, kg_bn], method="PBN", prob=0.5)
pbn = KGBN.load_network(pbn_string)

optimizer = KGBN.ParameterOptimizer(
    pbn,
    experiments=[
        {
            "perturbations": {"MYC": 1},
            "measurements": {"NRAS": 1.0},
        }
    ],
    config={
        "max_try": 1,
        "success_threshold": 0.1,
        "de_params": {"maxiter": 100, "popsize": 10, "polish": False, "workers": 1},
        "steady_state": {"method": "monte_carlo", "monte_carlo_params": {"n_runs": 2, "n_steps": 1000}},
        "seed": 9,
    },
)
result = optimizer.optimize(method='differential_evolution')
```

## Examples

Tutorials:

- **[BN_simulation.ipynb](./Examples/BN_simulation.ipynb)**: Basic Boolean network simulation
- **[PBN_simulation.ipynb](./Examples/PBN_simulation.ipynb)**: Probabilistic Boolean network simulation
- **[knowledge_graph.ipynb](./Examples/knowledge_graph.ipynb)**: Build from SIGNOR, extend, and merge into BN/PBN
- **[phenotype_score.ipynb](./Examples/phenotype_score.ipynb)**: Calculate phenotype scores using simulation results
- **[Optimization.ipynb](./Examples/Optimization.ipynb)**: Parameter optimization with experimental data
- **[BN_PBN_steady_state.ipynb](./Examples/BN_PBN_steady_state.ipynb)**: Steady state analysis
- **[BN_compression.ipynb](./Examples/BN_compression.ipynb)**: Model compression and simplification

Workflows:

- **[workflow_example.ipynb](./Examples/workflow_example.ipynb)**: Complete workflow from data to optimized model
- **[KGBN_workflow_toy.ipynb](./Examples/KGBN_workflow_toy.ipynb)**: A small toy workflow example

## Manuscript Reproduction

The manuscript reproduction notebook is kept in [manuscript/KGBN_manuscript.ipynb](manuscript/KGBN_manuscript.ipynb). 