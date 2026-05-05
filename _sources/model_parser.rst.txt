KGBN.model_parser
==================

The model_parser module provides functions for merging, converting, and extending Boolean and Probabilistic Boolean Networks.

.. automodule:: KGBN.model_parser
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
-------------

Merging Networks
~~~~~~~~~~~~~~~~

.. code-block:: python

   import KGBN

   # Load two networks
   network1 = KGBN.load_network_from_file("network1.txt")
   network2 = KGBN.load_network_from_file("network2.txt")

   # Merge into a Boolean-network definition string
   merged_bn_string = KGBN.merge_networks([network1, network2], method='Inhibitor Wins')
   merged_bn = KGBN.load_network(merged_bn_string, network_type='bn')

   # Or merge into a PBN definition string
   merged_pbn_string = KGBN.merge_networks([network1, network2], method='PBN', prob=0.9)
   merged_pbn = KGBN.load_network(merged_pbn_string, network_type='pbn')

Converting BN to PBN
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import KGBN

   # Load a Boolean Network
   bn = KGBN.load_network_from_file("network.txt")

   # Convert to PBN rules with equal probabilities for existing rules and self-loops
   pbn_string, nodes_to_optimize = KGBN.BN2PBN(bn, prob=0.5)
   pbn = KGBN.load_network(pbn_string, network_type='pbn')

Extending Networks
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import KGBN

   # Load original network and KG-derived network
   original_bn = KGBN.load_network_from_file("original.txt")
   kg_string, _ = KGBN.load_signor_network(gene_list=['TP53', 'MYC'])
   kg_network = KGBN.load_network_from_string(kg_string)

   # Extend original network with KG information
   extended_pbn_string = KGBN.extend_networks(
       original_bn, 
       kg_network, 
       nodes_to_extend=['GENE1', 'GENE2'],
       prob=0.5, # probability of the rules from the KG
       descriptive=True
   )
   extended_pbn = KGBN.load_network(extended_pbn_string, network_type='pbn')
