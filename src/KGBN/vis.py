
import igraph as ig
import re
import numpy as np
import json
import tempfile
import os
from pyvis.network import Network
import matplotlib.pyplot as plt
import networkx as nx

def is_entire_rule_negated(rule):
    """
    Check if the entire rule is a single negated expression like !(A | B).
    
    Returns True for rules like:
        - !(A | B)
        - !(A & B & C)
        - ! (A | B)
    
    Returns False for rules like:
        - !(A) | B  (negation is only partial)
        - A & !B    (negation is on individual variable)
        - !A        (no parentheses)
    """
    rule = rule.strip()
    
    # Check for !( or ! (
    if rule.startswith("!("):
        start_idx = 1
    elif rule.startswith("! ("):
        start_idx = 2
    else:
        return False
    
    # Check if opening paren's matching close paren is at the end
    depth = 0
    for i in range(start_idx, len(rule)):
        if rule[i] == '(':
            depth += 1
        elif rule[i] == ')':
            depth -= 1
            if depth == 0:
                # Found the matching closing paren
                return i == len(rule) - 1
    return False


def read_logic_rules(source):
    """
    Reads logic rules from a file path or from a string containing rules.

    Args:
        source (str): Path to the file or the string containing logic rules.

    Returns:
        dict: Mapping from variable names to their logic rules.
    """
    logic_rules = {}
    try:
        # Try to open as a file
        with open(source, 'r') as f:
            lines = f.readlines()
    except (OSError, TypeError):
        # If not a file, treat as string
        lines = source.splitlines()

    for line in lines:
        if line.startswith('#'):
            continue
        parts = line.strip().split('=')
        if len(parts) == 2:
            logic_rules[parts[0].strip()] = parts[1].strip()
    return logic_rules


def extract_logic_rules_from_network(network):
    """
    Extract logic rules from a BooleanNetwork or PBN object using stored equations.
    
    Args:
        network: BooleanNetwork or ProbabilisticBN object
        
    Returns:
        dict: Mapping from variable names to their logic rules (for BN) or list of rules (for PBN)
        dict: Mapping from edges to probabilities (for PBN) or empty dict (for BN)
    """
    logic_rules = {}
    edge_probabilities = {}
    
    if hasattr(network, 'equations') and network.equations:
        if hasattr(network, 'cij') and hasattr(network, 'gene_functions'):
            for equation in network.equations:
                if '=' in equation:
                    parts = equation.strip().split('=', 1)
                    if len(parts) == 2:
                        node_name = parts[0].strip()
                        rule = parts[1].strip()
                        
                        # For PBN, collect all rules for each node
                        if node_name not in logic_rules:
                            logic_rules[node_name] = []
                        logic_rules[node_name].append(rule)
            
            # Extract probability information for PBN
            for node_name, node_idx in network.nodeDict.items():
                if node_idx < len(network.nf):
                    num_funcs = network.nf[node_idx]
                    if num_funcs > 1:
                        # Multiple functions - extract probabilities
                        for func_offset in range(num_funcs):
                            if func_offset < len(network.cij[node_idx]):
                                probability = network.cij[node_idx, func_offset]
                                if probability > 0:
                                    # Map probability to specific rule
                                    if node_name in logic_rules and func_offset < len(logic_rules[node_name]):
                                        rule = logic_rules[node_name][func_offset]
                                        edge_probabilities[(node_name, rule)] = probability
                    else:
                        # Single function - probability is 1.0
                        if node_name in logic_rules and len(logic_rules[node_name]) > 0:
                            rule = logic_rules[node_name][0]
                            edge_probabilities[(node_name, rule)] = 1.0
        else:
            # This is a BN - use stored equations directly
            for equation in network.equations:
                if '=' in equation:
                    parts = equation.strip().split('=', 1)
                    if len(parts) == 2:
                        node_name = parts[0].strip()
                        rule = parts[1].strip()
                        logic_rules[node_name] = rule
    else:
        print("No logic rules provided.")
    
    return logic_rules, edge_probabilities


def build_igraph_pbn(logic_rules, edge_probabilities):
    """Build igraph for PBN with multiple rules per node."""
    g = ig.Graph(directed=True)
    node_names = set()
    
    # Collect all node names
    for node_name, rules in logic_rules.items():
        node_names.add(node_name)
        if isinstance(rules, list):
            for rule in rules:
                node_names.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
        else:
            node_names.update(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
    
    node_names = list(node_names)
    g.add_vertices(node_names)
    
    # Add edges for each rule
    for node_name, rules in logic_rules.items():
        if isinstance(rules, list):
            # PBN with multiple rules
            for rule in rules:
                inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
                prob = edge_probabilities.get((node_name, rule), 1.0)
                for input_node in inputs:
                    g.add_edge(input_node, node_name, label=rule, probability=prob)
        else:
            # BN with single rule
            inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
            for input_node in inputs:
                g.add_edge(input_node, node_name, label=rules, probability=1.0)
    
    return g


def build_igraph(logic_rules):
    """Build igraph for BN with single rule per node."""
    g = ig.Graph(directed=True)
    node_names = set(logic_rules.keys())

    # Also collect all appearing variables
    for rule in logic_rules.values():
        node_names.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))

    node_names = list(node_names)
    g.add_vertices(node_names)

    # Add directed edges
    for target, rule in logic_rules.items():
        inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
        for input_node in inputs:
            g.add_edge(input_node, target, label=rule)
    return g


def create_matplotlib_visualization(logic_rules, removed_nodes=None, removed_edges=None, 
                                  measured_nodes=None, perturbed_nodes=None, color_node=None):
    """
    Create a matplotlib-based visualization for Jupyter notebooks.
    
    Args:
        color_node (str): If provided, all nodes will use this color (overrides default coloring)
    """
    removed_nodes = removed_nodes or set()
    removed_edges = removed_edges or set()
    measured_nodes = measured_nodes or set()
    perturbed_nodes = perturbed_nodes or set()
    
    # Create networkx graph
    G = nx.DiGraph()
    
    # Add all nodes mentioned in rules
    all_nodes = set()
    if isinstance(list(logic_rules.values())[0], list):
        # PBN case
        for node_name, rules in logic_rules.items():
            all_nodes.add(node_name)
            for rule in rules:
                all_nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
    else:
        # BN case
        all_nodes = set(logic_rules.keys())
        for rule in logic_rules.values():
            all_nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
    
    # Add nodes with attributes
    for node in all_nodes:
        G.add_node(node, removed=node in removed_nodes)
    
    # Add edges
    for target, rules in logic_rules.items():
        if isinstance(rules, list):
            # PBN case
            for rule in rules:
                inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
                for input_node in inputs:
                    if input_node != target:
                        G.add_edge(input_node, target, rule=rule, removed=(input_node, target) in removed_edges)
        else:
            # BN case
            inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
            for input_node in inputs:
                if input_node != target:
                    G.add_edge(input_node, target, rule=rules, removed=(input_node, target) in removed_edges)
    
    # Create visualization
    plt.figure(figsize=(12, 8))
    
    # Use spring layout for better visualization
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # First, separate nodes into removed and normal nodes
    removed_node_list = [n for n in G.nodes() if n in removed_nodes]
    normal_nodes = [n for n in G.nodes() if n not in removed_nodes]
    
    # Categorize nodes based on provided sets and connectivity
    perturbed_node_list = [n for n in normal_nodes if n in perturbed_nodes]
    measured_node_list = [n for n in normal_nodes if n in measured_nodes and n not in perturbed_nodes]
    input_nodes = []
    output_nodes = []
    intermediate_nodes = []
    
    for node in normal_nodes:
        # Skip if already categorized as perturbed or measured
        if node in perturbed_nodes or node in measured_nodes:
            continue
            
        # Check if it's an input node (self-referential rule or no incoming edges)
        is_input = False
        
        if isinstance(logic_rules, dict):
            if isinstance(logic_rules.get(node), str):
                rule = logic_rules[node]
                if rule.strip() == node:
                    is_input = True
        
        # Also consider nodes with no incoming edges as inputs
        if G.in_degree(node) == 0:
            is_input = True
        
        # Check if node is output
        is_output = False
        if not measured_nodes:
            if G.out_degree(node) == 0 and not is_input:
                is_output = True
                

        # Categorize
        if is_input:
            input_nodes.append(node)
        elif is_output:
            output_nodes.append(node)
        else:
            intermediate_nodes.append(node)
    
    # Draw different node types with proper color scheme
    # If a uniform color is provided, use it for all nodes
    if color_node:
        all_normal_nodes = input_nodes + output_nodes + intermediate_nodes + perturbed_node_list + measured_node_list
        if all_normal_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=all_normal_nodes, node_color=color_node, 
                                 node_size=500, alpha=0.8, edgecolors='black')
        if removed_node_list:
            nx.draw_networkx_nodes(G, pos, nodelist=removed_node_list, node_color='lightgrey', 
                                 node_size=500, alpha=0.6, label='Removed', edgecolors='black')
    else:
        if input_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=input_nodes, node_color='lightgreen', 
                                 node_size=500, alpha=0.8, label='Input', edgecolors='black')
        if output_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=output_nodes, node_color='yellow', 
                                 node_size=500, alpha=0.8, label='Output', edgecolors='black')
        if intermediate_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=intermediate_nodes, node_color='lightblue', 
                                 node_size=500, alpha=0.8, label='Intermediate', edgecolors='black')
        if removed_node_list:
            nx.draw_networkx_nodes(G, pos, nodelist=removed_node_list, node_color='lightgrey', 
                                 node_size=500, alpha=0.6, label='Removed', edgecolors='black')
        if perturbed_node_list:
            nx.draw_networkx_nodes(G, pos, nodelist=perturbed_node_list, node_color='red', 
                                 node_size=500, alpha=0.8, label='Perturbed', edgecolors='black')
        if measured_node_list:
            nx.draw_networkx_nodes(G, pos, nodelist=measured_node_list, node_color='orange', 
                                 node_size=500, alpha=0.8, label='Measured', edgecolors='black')
    # Draw edges
    normal_edges = [(u, v) for u, v in G.edges() if (u, v) not in removed_edges]
    removed_edge_list = [(u, v) for u, v in G.edges() if (u, v) in removed_edges]
    
    if normal_edges:
        nx.draw_networkx_edges(G, pos, edgelist=normal_edges, edge_color='black', 
                             arrows=True, arrowsize=20, alpha=0.6)
    if removed_edge_list:
        nx.draw_networkx_edges(G, pos, edgelist=removed_edge_list, edge_color='gray', 
                             arrows=True, arrowsize=20, alpha=0.3, style='dashed')
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
    
    plt.title("Boolean Network Visualization", fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    
    # Add legend if there are different node types (not when uniform color is used)
    if not color_node and (perturbed_node_list or measured_node_list or input_nodes or intermediate_nodes or output_nodes or removed_node_list):
        plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    
    plt.show()
    return None


def vis_network(source, output_html="network_graph.html", interactive=False,
                removed_nodes=None, removed_edges=None, measured_nodes=None, perturbed_nodes=None,
                color_node=None, color_edge=None, physics=True,
                figsize=(13, 8), title="Boolean Network", node_groups=None,
                layout='hierarchical', layout_kwargs=None,
                activation_color='#D73027', inhibition_color='#4575B4',
                node_size=1500, font_size=9, ax=None, return_fig=False, seed=42):
    """
    Visualize the Boolean network.

    When ``interactive=False`` (default), produces a publication-quality
    matplotlib figure with a configurable layout and clearly distinguished
    activation / inhibition edges.

    When ``interactive=True``, saves an interactive HTML file via PyVis.

    Args:
        source: Logic rules — dict, file path (str), network string (str),
                BooleanNetwork, or ProbabilisticBN object.
        output_html (str): Output HTML file name (interactive mode only).
        interactive (bool): If ``True``, generate an interactive HTML file.
        removed_nodes (set): Nodes shown as removed (light grey).
        removed_edges (set): Edges shown as removed (interactive mode only).
        measured_nodes (set): Measured / readout nodes (orange).
        perturbed_nodes (set): Perturbed nodes (red).
        color_node (str): Uniform node colour override (interactive mode only).
        color_edge (str): Uniform edge colour override (interactive mode only).
        physics (bool): Enable physics simulation (interactive mode only).

        figsize (tuple): Figure size ``(width, height)`` in inches
                         (static mode only).
        title (str): Figure title (static mode only).
        node_groups (dict): Custom colour groups that override
                            connectivity-based defaults.  Format::

                                {group_name: (node_set, color)}

                            Example::

                                {'Drug targets': ({'BCL2', 'SYK'}, '#E1F2D0'),
                                 'Mutations':    ({'NRAS'},         '#D783FF')}

        layout (str): Node placement algorithm (static mode only).  Options:

                      * ``'hierarchical'`` *(default)* — layered top-down.
                      * ``'spring'`` — force-directed (dense/cyclic nets).
                      * ``'kamada_kawai'`` — stress-minimisation.
                      * ``'circular'`` — evenly spaced on a circle.
                      * ``'shell'`` — concentric shells.
                      * ``'spectral'`` — eigenvector-based.
                      * ``'dot'`` — Graphviz ``dot`` (requires
                        ``pygraphviz`` or ``pydot``).

        layout_kwargs (dict): Extra keyword arguments for the layout
                              algorithm.  For ``'hierarchical'``, ``x_gap``
                              (default 3.0) and ``y_gap`` (default 2.5)
                              control inter-node spacing.
        activation_color (str): Edge colour for activation (static mode only).
        inhibition_color (str): Edge colour for inhibition (static mode only).
        node_size (int): Node marker size in matplotlib units (static mode only).
        font_size (int): Label font size in points (static mode only).
        ax: Existing :class:`matplotlib.axes.Axes` to draw on; a new figure
            is created when *None* (static mode only).
        return_fig (bool): If ``True``, return ``(fig, ax)`` instead of
                           calling ``plt.show()`` (static mode only).
        seed (int): Random seed for reproducible layouts (static mode only).

    Returns:
        ``None``, or ``(fig, ax)`` when *return_fig* is ``True`` and
        ``interactive=False``.
    """
    # Handle different input types
    if hasattr(source, 'nodeDict'):  # BooleanNetwork or ProbabilisticBN
        logic_rules, edge_probabilities = extract_logic_rules_from_network(source)
    elif isinstance(source, dict):
        logic_rules = source
        edge_probabilities = {}
    else:  # String or file path
        logic_rules = read_logic_rules(source)
        edge_probabilities = {}

    # Check if logic_rules is empty
    if not logic_rules:
        print("No logic rules provided.")
        return

    # Set removed nodes and edges to empty sets if not provided
    removed_nodes = removed_nodes or set()
    removed_edges = removed_edges or set()
    measured_nodes = measured_nodes or set()
    perturbed_nodes = perturbed_nodes or set()

    # For non-interactive mode, use publication-quality matplotlib rendering
    if not interactive:
        return _vis_network_static(
            source,
            figsize=figsize, title=title,
            removed_nodes=removed_nodes, measured_nodes=measured_nodes,
            perturbed_nodes=perturbed_nodes, node_groups=node_groups,
            layout=layout, layout_kwargs=layout_kwargs,
            activation_color=activation_color, inhibition_color=inhibition_color,
            node_size=node_size, font_size=font_size,
            ax=ax, return_fig=return_fig, seed=seed,
        )

    # For interactive mode, use pyvis
    # Check if this is PBN (multiple rules per node)
    is_pbn = isinstance(list(logic_rules.values())[0], list) if logic_rules else False
    
    if is_pbn:
        # Build graph for PBN
        g = build_igraph_pbn(logic_rules, edge_probabilities)
    else:
        # Build graph for BN
        g = build_igraph(logic_rules)
    
    # Create pyvis network with proper configuration
    net = Network(
        height='1000px', 
        width='1200px', 
        bgcolor='#ffffff',
        font_color='black',
        directed=True
    )
    
    # Set network options
    if physics:
        net.set_options("""
        var options = {
            "configure": {
                "enabled": false
            },
            "edges": {
                "color": {
                    "inherit": true
                },
                "smooth": {
                    "enabled": true,
                    "type": "dynamic"
                }
            },
            "interaction": {
                "dragNodes": true,
                "hideEdgesOnDrag": false,
                "hideNodesOnDrag": false
            },
            "physics": {
                "enabled": true,
                "stabilization": {
                    "enabled": true,
                    "fit": true,
                    "iterations": 1000,
                    "onlyDynamicEdges": false,
                    "updateInterval": 50
                }
            }
        }
        """)
    else:
        net.set_options("""
        var options = {
            "configure": {
                "enabled": false
            },
            "edges": {
                "color": {
                    "inherit": true
                },
                "smooth": {
                    "enabled": true,
                    "type": "dynamic"
                }
            },
            "interaction": {
                "dragNodes": true,
                "hideEdgesOnDrag": false,
                "hideNodesOnDrag": false
            },
            "physics": {
                "enabled": false
            }
        }
        """)

    # Add nodes to the PyVis network
    nodes_added = []
    for v in g.vs:
        node_name = v["name"]
        upstream_count = len([pred for pred in g.predecessors(v.index) if pred != v.index])
        downstream_count = len([succ for succ in g.successors(v.index) if succ != v.index])
        
        # If uniform color is provided, use it for all non-removed nodes
        if color_node and node_name not in removed_nodes:
            node_color = color_node
            font_color = "black"
            nodes_added.append(f"{node_name}")
        # Set color based on removal status and node type
        elif node_name in removed_nodes:
            node_color = "lightgrey"
            font_color = "grey"
            nodes_added.append(f"{node_name} (removed)")
        elif node_name in perturbed_nodes:
            node_color = "red"
            font_color = "white"
            nodes_added.append(f"{node_name} (perturbed)")
        elif node_name in measured_nodes:
            node_color = "orange"
            font_color = "black"
            nodes_added.append(f"{node_name} (measured)")
        else:
            # Check if node is input based on logic rules or connectivity
            is_input = False
            
            # self-referential rule or no incoming edges
            if isinstance(logic_rules.get(node_name), str):
                rule = logic_rules[node_name]
                if rule.strip() == node_name:
                    is_input = True
            elif isinstance(logic_rules.get(node_name), list):
                # For PBN, check if all rules are self-referential
                rules = logic_rules[node_name]
                if all(rule.strip() == node_name for rule in rules):
                    is_input = True
            
            # Also consider nodes with no incoming edges as inputs
            if upstream_count == 0:
                is_input = True
            
            # Check if node is output
            is_output = False
            if not measured_nodes:
                if downstream_count == 0 and not is_input:
                    is_output = True
            
            # Set colors based on node type
            if is_input:
                node_color = "lightgreen"  # Light green for inputs
                font_color = "black"
                nodes_added.append(f"{node_name} (input)")
            elif is_output:
                node_color = "yellow"  # Yellow for outputs
                font_color = "black"
                nodes_added.append(f"{node_name} (output)")
            else:
                node_color = "lightblue"  # Light blue for intermediate nodes
                font_color = "black"
                nodes_added.append(f"{node_name} (intermediate)")

        net.add_node(
            v.index, 
            label=node_name, 
            title=node_name, 
            color=node_color, 
            size=40,
            font={'size': 20, 'color': font_color},
            shape='dot'
        )

    # Add edges to the PyVis network
    edges_added = []
    for e in g.es:
        src, tgt = e.tuple
        src_name = g.vs[src]["name"]
        tgt_name = g.vs[tgt]["name"]
        
        rule_label = e["label"]
        
        # Check if edge was removed
        edge_removed = (src_name, tgt_name) in removed_edges
        
        # Get probability if available
        edge_prob = e["probability"] if "probability" in e.attributes() else 1.0
        
        # Create hover title with rule and probability
        if is_pbn and edge_prob < 1.0:
            edge_title = f"{tgt_name} = {rule_label}, p={edge_prob:.3f}"
        else:
            edge_title = f"{tgt_name} = {rule_label}"
        
        # Set edge color and style based on rule type and removal status
        if edge_removed:
            edge_color = "lightgrey"
            edge_width = 1
            edge_alpha = 0.3
            edges_added.append(f"{src_name}->{tgt_name} (removed)")
        elif color_edge:
            # Use uniform edge color if provided
            edge_color = color_edge
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name}")
        elif is_entire_rule_negated(rule_label) and src_name in rule_label:
            # Entire rule is negated like !(A | B), all sources are direct inhibitors
            edge_color = "blue"
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name} (direct inhibitory)")
        elif f"!{src_name}" in rule_label or f"! {src_name}" in rule_label:
            # Direct negation of specific variable like !A or ! A
            edge_color = "blue"
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name} (direct inhibitory)")
        elif rule_label.startswith("!(") and src_name in rule_label:
            # Partial negation like !(A) | B, source is inside negated part
            edge_color = "grey"
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name} (inhibitory)")
        elif rule_label.startswith("! (") and src_name in rule_label:
            edge_color = "grey" 
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name} (inhibitory)")
        elif "! (" in rule_label and src_name in rule_label:
            edge_color = "grey"
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name} (inhibitory)")
        else:
            edge_color = "red"
            edge_width = 2
            edge_alpha = edge_prob if is_pbn else 1.0
            edges_added.append(f"{src_name}->{tgt_name} (activating)")
        
        # For PBN, adjust edge width and transparency based on probability
        if is_pbn and not edge_removed:
            edge_width = max(1, int(edge_prob * 4))  # Scale width by probability
            # Convert alpha to hex color for transparency
            alpha_hex = format(int(edge_alpha * 255), '02x')
            if edge_color == "red":
                edge_color = f"#ff0000{alpha_hex}"
            elif edge_color == "blue":
                edge_color = f"#0000ff{alpha_hex}"
            elif edge_color == "grey":
                edge_color = f"#808080{alpha_hex}"

        net.add_edge(
            src, 
            tgt, 
            title=edge_title, 
            color=edge_color, 
            width=edge_width
        )

    # Add legend for interactive mode
    legend_x = -1500  # Position legend to the left
    legend_y_start = -800
    legend_spacing = 120
    legend_idx = 0
    
    # Add node legend (skip if uniform color is used)
    if not color_node:
        legend_items = [
            ("Input", "lightgreen", "black"),
            ("Output", "yellow", "black"), 
            ("Intermediate", "lightblue", "black"),
            ("Measured", "orange", "black"),
            ("Perturbed", "red", "white"),
            ("Removed", "lightgrey", "grey")
        ]
        
        for label, legend_color, legend_font_color in legend_items:
            net.add_node(
                f"legend_{legend_idx}",
                label=label,
                title=f"{label} nodes",
                color=legend_color,
                size=30,
                shape='dot',
                x=legend_x,
                y=legend_y_start + legend_idx * legend_spacing,
                physics=False,  # Keep legend nodes fixed
                font={'size': 14, 'color': legend_font_color}
            )
            legend_idx += 1
    
    # Add edge legend using colored box with arrow symbol
    edge_legend_y_start = legend_y_start + legend_idx * legend_spacing + 60
    
    if not color_edge:
        # Create edge legend for activation/inhibition
        edge_legend_items = [
            ("→ Activation", "red"),
            ("→ Inhibition", "blue"),
        ]
    else:
        # Show the custom edge color in legend
        edge_legend_items = [
            ("→ Edge", color_edge),
        ]
    
    for edge_label, edge_legend_color in edge_legend_items:
        net.add_node(
            f"legend_edge_{legend_idx}",
            label=edge_label,
            title=f"{edge_label.replace('→ ', '')} edge",
            color=edge_legend_color,
            size=25,
            shape='box',
            x=legend_x,
            y=edge_legend_y_start,
            physics=False,
            font={'size': 14, 'color': 'white'}
        )
        edge_legend_y_start += legend_spacing
        legend_idx += 1

    # Save the network
    net.save_graph(output_html)
    print(f"Network visualization saved to {output_html}")


def vis_compression(original_network, compressed_network, 
                             compression_info, output_html="compression_comparison.html",
                             interactive=False):
    """
    Visualize the original network with removed/collapsed nodes highlighted.
    
    Args:
        original_network: Original BooleanNetwork or ProbabilisticBN
        compressed_network: Compressed network (not used for visualization)
        compression_info: Dictionary with compression information
        output_html (str): Output HTML file name
        interactive (bool): If True, return network visualization in interactive html file
    """
    # Collect all removed nodes from compression info
    removed_nodes = set()
    
    # Add explicitly removed nodes
    removed_nodes.update(compression_info.get('removed_non_observable', set()))
    removed_nodes.update(compression_info.get('removed_non_controllable', set()))
    
    # Add all removed nodes from the general tracking
    removed_nodes.update(compression_info.get('removed_nodes', set()))
    
    # Add intermediate nodes from collapsed paths
    for path in compression_info.get('collapsed_paths', []):
        if len(path) > 2:
            # All intermediate nodes (exclude first and last)
            removed_nodes.update(path[1:-1])
    
    # Get removed edges
    removed_edges = compression_info.get('removed_edges', set())
    
    measured_nodes = compression_info.get('measured_nodes', set())
    perturbed_nodes = compression_info.get('perturbed_nodes', set())
    
    # Use the original network for visualization with removed nodes marked
    return vis_network(
        original_network, 
        output_html, 
        interactive, 
        removed_nodes=removed_nodes,
        removed_edges=removed_edges,
        measured_nodes=measured_nodes,
        perturbed_nodes=perturbed_nodes
    )


def vis_extension(original_network, extended_network, output_html="network+KG.html", interactive=True,
                  color_node='#AED6F1', color_edge=None, extension_color_node='#FBE1BE',
                  extension_color_edge=None, physics=True,
                  figsize=(22, 9), title=None, node_groups=None,
                  layout='hierarchical', layout_kwargs=None, edge_scores=None,
                  activation_color='#D73027', inhibition_color='#4575B4',
                  node_size=1500, font_size=9, return_fig=False, seed=9):
    """
    Visualize the network with KG (Knowledge Graph) extension, highlighting
    new nodes and edges.

    When ``interactive=True`` (default), saves an interactive HTML file via
    PyVis.  When ``interactive=False``, produces a side-by-side
    publication-quality matplotlib figure comparing the original and extended
    networks.

    Args:
        original_network: Original BooleanNetwork or ProbabilisticBN.
        extended_network: Extended network with additional nodes / edges.
        output_html (str): Output HTML file name (interactive mode only).
        interactive (bool): If ``True``, generate an interactive HTML file.
        color_node (str): Fill colour for original-network nodes.
        color_edge (str): Uniform edge colour override (interactive mode only).
        extension_color_node (str): Fill colour for nodes added by the KG.
        extension_color_edge (str): Edge colour override for KG-added edges
                                    (interactive mode only).
        physics (bool): Enable physics simulation (interactive mode only).

        figsize (tuple): Overall figure size ``(width, height)`` in inches
                         (static mode only).
        title (str): Super-title above both panels (static mode only).
        node_groups (dict): Optional extra colour groups applied to both
                            panels, overriding *color_node* /
                            *extension_color_node* for the listed nodes.
                            Format::

                                {group_name: (node_set, color)}

                            Example::

                                {'Drug targets': ({'BCL2', 'SYK'}, '#FBE1BE'),
                                 'Mutations':    ({'NRAS'},         '#D783FF')}

        layout (str): Node placement algorithm shared by both panels
                      (static mode only).  See :func:`vis_network` for
                      available options.
        layout_kwargs (dict): Extra keyword arguments for the layout
                              algorithm (static mode only).
        edge_scores (dict): Optional ``{(src, tgt): score}`` mapping used to
                            scale new-edge line widths by evidence score
                            (static mode only).  Build this from the
                            ``all_relations`` list returned by
                            :func:`load_signor_network`::

                                _, relations = load_signor_network(genes, ...)
                                edge_scores = {
                                    (src, tgt): score
                                    for src, tgt, _, score in relations
                                    if score is not None
                                }

                            Width is mapped as ``0.5 + score × 4``.
        activation_color (str): Colour for activation edges (static mode only).
        inhibition_color (str): Colour for inhibition edges (static mode only).
        node_size (int): Node marker size (static mode only).
        font_size (int): Label font size in points (static mode only).
        return_fig (bool): If ``True``, return ``(fig, axes)`` instead of
                           calling ``plt.show()`` (static mode only).
        seed (int): Random seed for reproducible layouts (static mode only).

    Returns:
        ``None``, or ``(fig, axes)`` when *return_fig* is ``True`` and
        ``interactive=False``.
    """
    # Extract logic rules from both networks
    original_rules, _ = extract_logic_rules_from_network(original_network)
    extended_rules, extended_probabilities = extract_logic_rules_from_network(extended_network)
    
    # Helper function to extract nodes from rules
    def extract_nodes_from_rules(rules_dict):
        nodes = set()
        for node_name, rules in rules_dict.items():
            nodes.add(node_name)
            if isinstance(rules, list):
                # PBN case - rules is a list of rule strings
                for rule in rules:
                    nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
            else:
                # BN case - rules is a single string
                nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
        return nodes
    
    # Helper function to extract edges from rules
    def extract_edges_from_rules(rules_dict):
        edges = set()
        for target, rules in rules_dict.items():
            if isinstance(rules, list):
                # PBN case
                for rule in rules:
                    inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
                    for input_node in inputs:
                        edges.add((input_node, target))
            else:
                # BN case
                inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
                for input_node in inputs:
                    edges.add((input_node, target))
        return edges
    
    # Extract nodes and edges
    original_nodes = extract_nodes_from_rules(original_rules)
    extended_nodes = extract_nodes_from_rules(extended_rules)
    original_edges = extract_edges_from_rules(original_rules)
    extended_edges = extract_edges_from_rules(extended_rules)
    
    # Identify new nodes and edges
    new_nodes = extended_nodes - original_nodes
    new_edges = extended_edges - original_edges
    
    print(f"Extension comparison:")
    print(f"  Original nodes: {len(original_nodes)}")
    print(f"  Extended nodes: {len(extended_nodes)}")
    print(f"  New nodes: {len(new_nodes)} - {sorted(new_nodes)}")
    print(f"  Original edges: {len(original_edges)}")
    print(f"  Extended edges: {len(extended_edges)}")
    print(f"  New edges: {len(new_edges)}")
    
    # For non-interactive mode, use publication-quality side-by-side rendering
    if not interactive:
        return _vis_extension_static(
            original_network, extended_network,
            figsize=figsize, title=title,
            color_node=color_node, extension_color_node=extension_color_node,
            node_groups=node_groups,
            layout=layout, layout_kwargs=layout_kwargs,
            edge_scores=edge_scores,
            activation_color=activation_color, inhibition_color=inhibition_color,
            node_size=node_size, font_size=font_size,
            return_fig=return_fig, seed=seed,
        )
    
    # For interactive mode, use the extended network but highlight new elements
    # Check if this is PBN
    is_pbn = any(isinstance(rules, list) for rules in extended_rules.values()) if extended_rules else False
    
    if is_pbn:
        # Build graph for PBN
        g = build_igraph_pbn(extended_rules, extended_probabilities)
    else:
        # Build graph for BN
        g = build_igraph(extended_rules)
    
    # Create pyvis network
    net = Network(
        height='1000px', 
        width='1200px', 
        bgcolor='#ffffff',
        font_color='black',
        directed=True
    )
    
    # Set network options
    if physics:
        net.set_options("""
        var options = {
            "configure": {
                "enabled": false
            },
            "edges": {
                "color": {
                    "inherit": true
                },
                "smooth": {
                    "enabled": true,
                    "type": "dynamic"
                }
            },
            "interaction": {
                "dragNodes": true,
                "hideEdgesOnDrag": false,
                "hideNodesOnDrag": false
            },
            "physics": {
                "enabled": true,
                "stabilization": {
                    "enabled": true,
                    "fit": true,
                    "iterations": 1000,
                    "onlyDynamicEdges": false,
                    "updateInterval": 50
                }
            }
        }
        """)
    else:
        net.set_options("""
        var options = {
            "configure": {
                "enabled": false
            },
            "edges": {
                "color": {
                    "inherit": true
                },
                "smooth": {
                    "enabled": true,
                    "type": "dynamic"
                }
            },
            "interaction": {
                "dragNodes": true,
                "hideEdgesOnDrag": false,
                "hideNodesOnDrag": false
            },
            "physics": {
                "enabled": false
            }
        }
        """)

    # Add nodes to the PyVis network
    for v in g.vs:
        node_name = v["name"]
        
        # Check if this is a new node (from KG extension)
        is_new_node = node_name in new_nodes
        
        # Set node color based on whether it's from original network or KG extension
        if is_new_node:
            node_color = extension_color_node
            title = f"{node_name} (from KG)"
        elif color_node:
            node_color = color_node
            title = f"{node_name} (original)"
        else:
            # Default color based on node type
            node_color = "lightblue"
            title = f"{node_name} (original)"
        
        font_color = "black"

        net.add_node(
            v.index, 
            label=node_name, 
            title=title, 
            color=node_color,
            size=40,
            font={'size': 20, 'color': font_color},
            shape='dot'
        )

    # Add edges to the PyVis network
    for e in g.es:
        src, tgt = e.tuple
        src_name = g.vs[src]["name"]
        tgt_name = g.vs[tgt]["name"]
        
        rule_label = e["label"]
        
        # Check if this is a new edge
        is_new_edge = (src_name, tgt_name) in new_edges
        
        # Get probability if available
        edge_prob = e["probability"] if "probability" in e.attributes() else 1.0
        
        # Create hover title with rule and probability
        if is_pbn and edge_prob < 1.0:
            edge_title = f"{tgt_name} = {rule_label}, p={edge_prob:.3f}"
        else:
            edge_title = f"{tgt_name} = {rule_label}"
        
        if is_new_edge:
            edge_title += " (from KG)"
        
        # Determine if this is an inhibitory edge based on rule
        is_inhibitory = (
            (is_entire_rule_negated(rule_label) and src_name in rule_label) or
            f"!{src_name}" in rule_label or 
            f"! {src_name}" in rule_label or
            (rule_label.startswith("!(") and src_name in rule_label) or
            (rule_label.startswith("! (") and src_name in rule_label) or
            ("! (" in rule_label and src_name in rule_label)
        )
        
        # Set edge color and style
        if is_new_edge:
            # New edge from KG
            if extension_color_edge:
                edge_color = extension_color_edge
            else:
                # Fall back to regulation type
                edge_color = "blue" if is_inhibitory else "red"
            edge_width = 3
        else:
            # Original network edge
            if color_edge:
                edge_color = color_edge
            else:
                # Fall back to regulation type
                edge_color = "blue" if is_inhibitory else "red"
            edge_width = 2
        
        # For PBN, adjust edge width and transparency based on probability
        if is_pbn and not is_new_edge and edge_prob < 1.0:
            edge_width = max(1, int(edge_prob * 4))  # Scale width by probability
            # Convert alpha to hex color for transparency
            alpha_hex = format(int(edge_prob * 255), '02x')
            if edge_color == "red":
                edge_color = f"#ff0000{alpha_hex}"
            elif edge_color == "blue":
                edge_color = f"#0000ff{alpha_hex}"
            elif edge_color == "grey":
                edge_color = f"#808080{alpha_hex}"

        net.add_edge(
            src, 
            tgt, 
            title=edge_title, 
            color=edge_color, 
            width=edge_width
        )

    # Add legend
    legend_x = -1500  
    legend_y_start = -800
    legend_spacing = 120
    legend_idx = 0
    
    # Build node legend items based on what colors are being used
    node_legend_items = []
    if color_node:
        node_legend_items.append(("Original Network", color_node, "black"))
    node_legend_items.append(("KG Extension", extension_color_node, "black"))
    
    for label, legend_color, font_color in node_legend_items:
        net.add_node(
            f"legend_{legend_idx}",
            label=label,
            title=f"{label}",
            color=legend_color,
            size=30,
            shape='dot',
            x=legend_x,
            y=legend_y_start + legend_idx * legend_spacing,
            physics=False,
            font={'size': 14, 'color': font_color}
        )
        legend_idx += 1
    
    # Add edge legend using colored box with arrow symbol
    edge_legend_y_start = legend_y_start + legend_idx * legend_spacing + 60
    
    # Build edge legend based on what colors are being used
    edge_legend_items = []
    
    # Original network edges legend
    if color_edge:
        edge_legend_items.append(("→ Original Edge", color_edge))
    else:
        # Original edges use regulation-based colors
        edge_legend_items.append(("→ Activation", "red"))
        edge_legend_items.append(("→ Inhibition", "blue"))
    
    # KG extension edges legend
    if extension_color_edge:
        edge_legend_items.append(("→ KG Edge", extension_color_edge))
    elif color_edge:
        # KG edges use regulation-based colors (only add if not already shown above)
        edge_legend_items.append(("→ KG Activation", "red"))
        edge_legend_items.append(("→ KG Inhibition", "blue"))
    
    for edge_label, edge_legend_color in edge_legend_items:
        net.add_node(
            f"legend_edge_{legend_idx}",
            label=edge_label,
            title=f"{edge_label.replace('→ ', '')} edge",
            color=edge_legend_color,
            size=25,
            shape='box',
            x=legend_x,
            y=edge_legend_y_start,
            physics=False,
            font={'size': 14, 'color': 'white'}
        )
        edge_legend_y_start += legend_spacing
        legend_idx += 1

    # Save the network
    net.save_graph(output_html)
    print(f"Extension visualization saved to {output_html}")


def create_matplotlib_extension_visualization(logic_rules, new_nodes, new_edges,
                                              color_node='lightblue', extension_color_node='#E1F2D0',
                                              extension_color_edge=None):
    """
    Create a matplotlib-based visualization for extension comparison.
    
    Args:
        logic_rules: Dictionary of logic rules
        new_nodes: Set of new node names from KG extension
        new_edges: Set of new edge tuples from KG extension
        color_node (str): Color for nodes from the original network
        extension_color_node (str): Color for new nodes from KG
        extension_color_edge (str): Color for new edges from KG. If None, uses default color
    """
    # Create networkx graph
    G = nx.DiGraph()
    
    # Add all nodes mentioned in rules
    all_nodes = set()
    for node_name, rules in logic_rules.items():
        all_nodes.add(node_name)
        if isinstance(rules, list):
            # PBN case
            for rule in rules:
                all_nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
        else:
            # BN case
            all_nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
    
    # Add nodes with attributes
    for node in all_nodes:
        G.add_node(node, is_new=node in new_nodes)
    
    # Add edges
    for target, rules in logic_rules.items():
        if isinstance(rules, list):
            # PBN case
            for rule in rules:
                inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
                for input_node in inputs:
                    if input_node != target:
                        G.add_edge(input_node, target, is_new=(input_node, target) in new_edges)
        else:
            # BN case
            inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rules))
            for input_node in inputs:
                if input_node != target:
                    G.add_edge(input_node, target, is_new=(input_node, target) in new_edges)
    
    # Create visualization
    plt.figure(figsize=(12, 8))
    
    # Use spring layout
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Separate new and existing nodes
    existing_nodes = [n for n in G.nodes() if n not in new_nodes]
    new_node_list = [n for n in G.nodes() if n in new_nodes]
    
    # Draw nodes with custom colors
    original_node_color = color_node if color_node else 'lightblue'
    if existing_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=existing_nodes, node_color=original_node_color, 
                             node_size=500, alpha=0.8, label='Original Network', edgecolors='black')
    if new_node_list:
        nx.draw_networkx_nodes(G, pos, nodelist=new_node_list, node_color=extension_color_node, 
                             node_size=500, alpha=0.8, label='KG Extension', edgecolors='black')
    
    # Draw edges
    existing_edges = [(u, v) for u, v in G.edges() if (u, v) not in new_edges]
    new_edge_list = [(u, v) for u, v in G.edges() if (u, v) in new_edges]
    
    if existing_edges:
        nx.draw_networkx_edges(G, pos, edgelist=existing_edges, edge_color='black', 
                             arrows=True, arrowsize=20, alpha=0.6)
    if new_edge_list:
        new_edge_color = extension_color_edge if extension_color_edge else 'orange'
        nx.draw_networkx_edges(G, pos, edgelist=new_edge_list, edge_color=new_edge_color, 
                             arrows=True, arrowsize=20, alpha=0.8, width=3)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
    
    plt.title("Network + KG Extension Visualization", fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()

    # Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))

    plt.show()
    return None


# ---------------------------------------------------------------------------
# Publication-quality figure functions
# ---------------------------------------------------------------------------

def _get_edge_type(src_name, rule):
    """Determine if an edge from src_name is inhibitory based on the logic rule."""
    if is_entire_rule_negated(rule):
        tokens = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
        if src_name in tokens:
            return 'inhibitory'
    if f'!{src_name}' in rule or f'! {src_name}' in rule:
        return 'inhibitory'
    return 'activating'


def _build_nx_graph_from_rules(logic_rules):
    """Build a networkx DiGraph from logic rules dictionary, with edge type attributes."""
    G = nx.DiGraph()

    all_nodes = set()
    for node_name, rules in logic_rules.items():
        all_nodes.add(node_name)
        rule_list = rules if isinstance(rules, list) else [rules]
        for rule in rule_list:
            all_nodes.update(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
    G.add_nodes_from(sorted(all_nodes))

    for target, rules in sorted(logic_rules.items()):
        rule_list = rules if isinstance(rules, list) else [rules]
        for rule in rule_list:
            inputs = set(re.findall(r'\b[A-Za-z0-9_]+\b', rule))
            for src in sorted(inputs):
                if src == target:
                    continue
                edge_type = _get_edge_type(src, rule)
                if not G.has_edge(src, target):
                    G.add_edge(src, target, type=edge_type)
                elif G[src][target]['type'] != edge_type:
                    G[src][target]['type'] = 'ambiguous'

    return G


def _dag_from_graph(G):
    """
    Return a DAG derived from G by removing self-loops and back-edges found
    via an iterative DFS.  Forward/cross edges are preserved.
    """
    H = nx.DiGraph()
    H.add_nodes_from(G.nodes())
    for u, v in G.edges():
        if u != v:
            H.add_edge(u, v)

    if nx.is_directed_acyclic_graph(H):
        return H

    # Iterative DFS: skip edges that would create a back-edge (ancestor → descendant)
    dag = nx.DiGraph()
    dag.add_nodes_from(H.nodes())
    visited = set()
    on_stack = set()

    for start in sorted(H.nodes()):
        if start in visited:
            continue
        stack = [(start, iter(sorted(H.successors(start))))]
        visited.add(start)
        on_stack.add(start)
        while stack:
            node, children = stack[-1]
            try:
                child = next(children)
                if child in on_stack:
                    pass  # back edge — omit
                else:
                    dag.add_edge(node, child)
                    if child not in visited:
                        visited.add(child)
                        on_stack.add(child)
                        stack.append((child, iter(sorted(H.successors(child)))))
            except StopIteration:
                on_stack.discard(node)
                stack.pop()

    return dag


def _compute_hierarchical_pos(G, x_gap=3.0, y_gap=2.5):
    """
    Compute layered (Sugiyama-style) positions for a directed graph.
    Cycles are handled by removing DFS back-edges before layer assignment,
    so no node is pushed to an unexpectedly deep layer.

    Returns:
        dict: {node: np.array([x, y])} positions
    """
    if not G.nodes():
        return {}

    dag = _dag_from_graph(G)

    # Assign layers via longest path from sources on the DAG
    topo = list(nx.lexicographical_topological_sort(dag, key=lambda n: str(n)))
    layers = {n: 0 for n in dag.nodes()}
    for node in topo:
        for succ in dag.successors(node):
            if layers[succ] <= layers[node]:
                layers[succ] = layers[node] + 1

    # Group nodes by layer
    layer_groups = {}
    for node, layer in layers.items():
        layer_groups.setdefault(layer, []).append(node)

    max_layer = max(layer_groups.keys()) if layer_groups else 0

    # Barycenter ordering within layers to reduce edge crossings (3 passes)
    for _ in range(3):
        for layer_idx in sorted(layer_groups.keys()):
            nodes = layer_groups[layer_idx]
            if layer_idx == 0:
                nodes.sort(key=lambda n: -dag.out_degree(n))
            else:
                prev_layer = layer_groups.get(layer_idx - 1, [])
                prev_x = {n: i for i, n in enumerate(prev_layer)}

                def bary_key(n, pxmap=prev_x, pl=prev_layer):
                    preds = [p for p in dag.predecessors(n) if p in pxmap]
                    return np.mean([pxmap[p] for p in preds]) if preds else len(pl) / 2

                nodes.sort(key=bary_key)
            layer_groups[layer_idx] = nodes

    # Assign coordinates: top layer = highest y
    pos = {}
    for layer_idx, nodes in layer_groups.items():
        n = len(nodes)
        y = (max_layer - layer_idx) * y_gap
        for i, node in enumerate(nodes):
            x = (i - (n - 1) / 2) * x_gap
            pos[node] = np.array([x, y])

    return pos


def _compute_layout(G, layout='hierarchical', layout_kwargs=None):
    """
    Compute node positions using the requested algorithm.

    Args:
        G: networkx DiGraph
        layout (str): One of ``'hierarchical'``, ``'spring'``,
                      ``'kamada_kawai'``, ``'circular'``, ``'shell'``,
                      ``'spectral'``, or ``'dot'`` (requires graphviz).
        layout_kwargs (dict): Extra keyword arguments forwarded to the
                              chosen algorithm.  For ``'hierarchical'``
                              the supported keys are ``x_gap`` and
                              ``y_gap`` (floats, default 3.0 / 2.5).

    Returns:
        dict: {node: position} mapping accepted by networkx draw functions.
    """
    kw = layout_kwargs or {}

    if layout == 'hierarchical':
        return _compute_hierarchical_pos(
            G,
            x_gap=kw.get('x_gap', 3.0),
            y_gap=kw.get('y_gap', 2.5),
        )
    elif layout == 'spring':
        defaults = dict(k=2.0, iterations=100, seed=42)
        defaults.update(kw)
        return nx.spring_layout(G, **defaults)
    elif layout == 'kamada_kawai':
        return nx.kamada_kawai_layout(G, **kw)
    elif layout == 'circular':
        return nx.circular_layout(G, **kw)
    elif layout == 'shell':
        return nx.shell_layout(G, **kw)
    elif layout == 'spectral':
        return nx.spectral_layout(G, **kw)
    elif layout == 'dot':
        try:
            return nx.drawing.nx_agraph.graphviz_layout(G, prog='dot', **kw)
        except Exception:
            try:
                return nx.drawing.nx_pydot.graphviz_layout(G, prog='dot', **kw)
            except Exception:
                print("graphviz unavailable — falling back to hierarchical layout.")
                return _compute_hierarchical_pos(G)
    else:
        raise ValueError(
            f"Unknown layout '{layout}'. Choose from: 'hierarchical', "
            "'spring', 'kamada_kawai', 'circular', 'shell', 'spectral', 'dot'."
        )


def _draw_network_on_ax(ax, G, pos, node_color_map,
                        activation_color, inhibition_color,
                        node_size, font_size,
                        new_node_set=None, new_edge_set=None,
                        gray_edge_set=None, removed_nodes=None,
                        edge_scores=None,
                        missing_node_set=None, missing_edge_set=None):
    """
    Internal helper: draw graph G on ax with hierarchical style.

    new_edge_set     — edges drawn thick and fully colored (KG additions).
    gray_edge_set    — edges drawn thin and gray (solid=activation, dashed=inhibition).
    missing_node_set — nodes present in original but absent from extended network;
                       drawn with a dashed border to signal removal.
    missing_edge_set — edges present in original but absent from extended network;
                       drawn as light gray dotted lines.
    edge_scores      — dict {(src, tgt): score} to scale new-edge widths.
    """
    new_node_set    = new_node_set    or set()
    new_edge_set    = new_edge_set    or set()
    gray_edge_set   = gray_edge_set   or set()
    removed_nodes   = removed_nodes   or set()
    edge_scores     = edge_scores     or {}
    missing_node_set = missing_node_set or set()
    missing_edge_set = missing_edge_set or set()

    # Categorise edges into buckets
    act_base, inh_base, act_new, inh_new, gray_act, gray_inh, miss_act, miss_inh, faded = \
        [], [], [], [], [], [], [], [], []
    for u, v, d in G.edges(data=True):
        if u == v:
            continue
        if u in removed_nodes or v in removed_nodes:
            faded.append((u, v))
            continue
        etype = d.get('type', 'activating')
        if (u, v) in new_edge_set:
            (inh_new  if etype == 'inhibitory' else act_new ).append((u, v))
        elif (u, v) in missing_edge_set:
            (miss_inh if etype == 'inhibitory' else miss_act).append((u, v))
        elif (u, v) in gray_edge_set:
            (gray_inh if etype == 'inhibitory' else gray_act).append((u, v))
        else:
            (inh_base if etype == 'inhibitory' else act_base).append((u, v))

    base_kw = dict(ax=ax, pos=pos, node_size=node_size, arrows=True, arrowsize=15)

    def _new_widths(edgelist, default=3.0):
        if not edge_scores:
            return default
        return [1.5 + edge_scores.get((u, v), 0.5) * 4.0 for u, v in edgelist]

    for edgelist, color, style, width, alpha in [
        # Gray original edges: activation=solid, inhibition=dashed
        (gray_act,  '#AAAAAA', 'solid',  1.0, 0.55),
        (gray_inh,  '#AAAAAA', 'dashed', 1.0, 0.55),
        # Missing edges: lighter dotted, activation=solid dotted, inhibition=dashed dotted
        (miss_act,  '#CCCCCC', 'dotted', 1.0, 0.5),
        (miss_inh,  '#CCCCCC', 'dotted', 1.0, 0.5),
        # Colored base edges
        (act_base,  activation_color, 'solid',  1.5, 0.75),
        (inh_base,  inhibition_color, 'dashed', 1.5, 0.75),
        (faded,     '#DDDDDD', 'dotted', 1.0, 0.3),
    ]:
        if edgelist:
            nx.draw_networkx_edges(G, edgelist=edgelist, edge_color=color,
                                   style=style, width=width, alpha=alpha, **base_kw)

    # New edges: drawn separately so width can be score-scaled
    for edgelist, color, style, alpha in [
        (act_new, activation_color, 'solid',  1.0),
        (inh_new, inhibition_color, 'dashed', 1.0),
    ]:
        if edgelist:
            nx.draw_networkx_edges(G, edgelist=edgelist, edge_color=color,
                                   style=style, width=_new_widths(edgelist),
                                   alpha=alpha, **base_kw)

    # Draw nodes — split into present and missing so we can apply dashed borders
    node_list = list(G.nodes())
    present_nodes  = [n for n in node_list if n not in missing_node_set]
    absent_nodes   = [n for n in node_list if n in  missing_node_set]

    def _node_style(nodelist):
        colors  = [node_color_map.get(n, '#AED6F1') for n in nodelist]
        borders = ['#FFC000' if n in new_node_set else '#333333' for n in nodelist]
        lws     = [3.0 if n in new_node_set else 1.5 for n in nodelist]
        return colors, borders, lws

    if present_nodes:
        colors, borders, lws = _node_style(present_nodes)
        nx.draw_networkx_nodes(G, pos, nodelist=present_nodes, node_color=colors,
                               node_size=node_size, edgecolors=borders,
                               linewidths=lws, ax=ax)

    if absent_nodes:
        colors, borders, lws = _node_style(absent_nodes)
        # Draw with lighter fill to signal absence; dashed border via PathCollection
        faded_colors = [c + '99' if isinstance(c, str) and c.startswith('#') and len(c) == 7
                        else c for c in colors]
        pc = nx.draw_networkx_nodes(G, pos, nodelist=absent_nodes,
                                    node_color=faded_colors,
                                    node_size=node_size, edgecolors='#999999',
                                    linewidths=1.5, ax=ax)
        if pc is not None:
            pc.set_linestyle('dashed')

    # Labels
    active_labels  = {n: n for n in node_list if n not in removed_nodes}
    removed_labels = {n: n for n in node_list if n in  removed_nodes}
    if active_labels:
        nx.draw_networkx_labels(G, pos, labels=active_labels,
                                font_size=font_size, font_weight='bold', ax=ax)
    if removed_labels:
        nx.draw_networkx_labels(G, pos, labels=removed_labels,
                                font_size=font_size, font_color='#888888', ax=ax)


def _vis_network_static(source, figsize=(13, 8), title="Boolean Network",
                        removed_nodes=None, measured_nodes=None, perturbed_nodes=None,
                        node_groups=None,
                        layout='hierarchical', layout_kwargs=None,
                        activation_color='#D73027', inhibition_color='#4575B4',
                        node_size=1500, font_size=9,
                        ax=None, return_fig=False, seed=42):
    """
    Internal helper: publication-quality static network figure via matplotlib.

    Called by :func:`vis_network` when ``interactive=False``.
    """
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import random

    random.seed(seed)
    np.random.seed(seed)

    # --- Parse source ---
    if hasattr(source, 'nodeDict'):
        logic_rules, _ = extract_logic_rules_from_network(source)
    elif isinstance(source, dict):
        logic_rules = source
    else:
        logic_rules = read_logic_rules(source)

    if not logic_rules:
        print("No logic rules provided.")
        return

    removed_nodes = removed_nodes or set()
    measured_nodes = measured_nodes or set()
    perturbed_nodes = perturbed_nodes or set()
    node_groups = node_groups or {}

    # --- Build graph & layout ---
    G = _build_nx_graph_from_rules(logic_rules)
    layout_kwargs = dict(layout_kwargs or {})
    if layout == 'spring' and 'seed' not in layout_kwargs:
        layout_kwargs['seed'] = seed
    pos = _compute_layout(G, layout=layout, layout_kwargs=layout_kwargs)

    # --- Node colours (priority: perturbed > measured > groups > connectivity) ---
    node_color_map = {}
    for node in G.nodes():
        has_pred = any(p != node for p in G.predecessors(node))
        has_succ = any(s != node for s in G.successors(node))
        if not has_pred:
            node_color_map[node] = '#A8D5A2'   # light green — input
        elif not has_succ:
            node_color_map[node] = '#FFE484'   # yellow — output
        else:
            node_color_map[node] = '#AED6F1'   # light blue — intermediate

    for _, (group_nodes, group_color) in node_groups.items():
        for node in group_nodes:
            if node in node_color_map:
                node_color_map[node] = group_color

    for node in removed_nodes:
        node_color_map[node] = '#DDDDDD'
    for node in measured_nodes:
        if node not in removed_nodes:
            node_color_map[node] = '#F0A500'
    for node in perturbed_nodes:
        if node not in removed_nodes:
            node_color_map[node] = '#E74C3C'

    # --- Axes ---
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # --- Draw ---
    _draw_network_on_ax(ax, G, pos, node_color_map,
                        activation_color, inhibition_color,
                        node_size, font_size,
                        removed_nodes=removed_nodes)

    # --- Legend ---
    handles = []
    if node_groups:
        for group_name, (_, gc) in node_groups.items():
            handles.append(Patch(facecolor=gc, edgecolor='#333333', label=group_name))
    else:
        handles += [
            Patch(facecolor='#A8D5A2', edgecolor='#333333', label='Input'),
            Patch(facecolor='#FFE484', edgecolor='#333333', label='Output'),
            Patch(facecolor='#AED6F1', edgecolor='#333333', label='Intermediate'),
        ]
    if measured_nodes:
        handles.append(Patch(facecolor='#F0A500', edgecolor='#333333', label='Measured'))
    if perturbed_nodes:
        handles.append(Patch(facecolor='#E74C3C', edgecolor='#333333', label='Perturbed'))
    if removed_nodes:
        handles.append(Patch(facecolor='#DDDDDD', edgecolor='#AAAAAA', label='Removed'))
    handles += [
        Line2D([0], [0], color=activation_color, lw=2, ls='-',  label='Activation'),
        Line2D([0], [0], color=inhibition_color,  lw=2, ls='--', label='Inhibition'),
    ]

    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1),
              borderaxespad=0, fontsize=8, framealpha=0.85)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.axis('off')
    plt.tight_layout()

    if return_fig:
        return fig, ax
    plt.show()
    return None


def _vis_extension_static(original_network, extended_network,
                          figsize=(22, 9), title=None,
                          color_node='#AED6F1',
                          extension_color_node='#FBE1BE',
                          node_groups=None,
                          layout='hierarchical', layout_kwargs=None,
                          edge_scores=None,
                          activation_color='#D73027',
                          inhibition_color='#4575B4',
                          node_size=1500, font_size=9,
                          return_fig=False,
                          seed=9):
    """
    Internal helper: side-by-side publication figure comparing original and
    KG-extended networks via matplotlib.

    Called by :func:`vis_extension` when ``interactive=False``.
    """
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import random
    random.seed(seed)
    np.random.seed(seed)

    # --- Extract rules ---
    original_rules, _ = extract_logic_rules_from_network(original_network)
    extended_rules, _ = extract_logic_rules_from_network(extended_network)
    node_groups = node_groups or {}

    # --- Build graphs ---
    G_orig = _build_nx_graph_from_rules(original_rules)
    G_ext = _build_nx_graph_from_rules(extended_rules)

    orig_nodes = set(G_orig.nodes())
    ext_nodes = set(G_ext.nodes())
    orig_edges = set(G_orig.edges())
    ext_edges = set(G_ext.edges())
    new_nodes     = ext_nodes  - orig_nodes   # added by KG
    new_edges     = ext_edges  - orig_edges   # added by KG
    missing_nodes = orig_nodes - ext_nodes    # in original but dropped from extended
    missing_edges = orig_edges - ext_edges    # in original but dropped from extended

    print("Extension comparison:")
    print(f"  Original : {len(orig_nodes)} nodes, {len(orig_edges)} edges")
    print(f"  Extended : {len(ext_nodes)} nodes, {len(ext_edges)} edges")
    print(f"  New nodes ({len(new_nodes)}): {sorted(new_nodes)}")
    print(f"  New edges: {len(new_edges)}")
    if missing_nodes:
        print(f"  Missing nodes ({len(missing_nodes)}): {sorted(missing_nodes)}")
    if missing_edges:
        print(f"  Missing edges: {len(missing_edges)}")

    # --- Shared layout ---
    # Build a union graph (all nodes from both networks) so that missing nodes
    # (in original but not extended) also get positions, keeping shared nodes
    # at the same (x, y) across both panels.
    G_union = nx.DiGraph()
    G_union.add_nodes_from(G_orig.nodes())
    G_union.add_nodes_from(G_ext.nodes())
    G_union.add_edges_from(G_orig.edges())
    G_union.add_edges_from(G_ext.edges())
    layout_kwargs = dict(layout_kwargs or {})
    if layout == 'spring' and 'seed' not in layout_kwargs:
        layout_kwargs['seed'] = seed
    pos_shared = _compute_layout(G_union, layout=layout, layout_kwargs=layout_kwargs)
    pos_orig = {n: pos_shared[n] for n in G_orig.nodes()}

    # --- Colour maps ---
    def make_color_map(G, is_ext):
        cmap = {}
        for node in G.nodes():
            if is_ext and node in new_nodes:
                cmap[node] = extension_color_node
            else:
                cmap[node] = color_node
        for _, (group_nodes, group_color) in node_groups.items():
            for node in group_nodes:
                if node in cmap:
                    cmap[node] = group_color
        return cmap

    orig_cmap = make_color_map(G_orig, is_ext=False)
    ext_cmap = make_color_map(G_ext, is_ext=True)

    # --- Edge scores (for weighted mode) ---
    edge_scores = edge_scores or {}

    # --- Figure ---
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Left panel: all edges shown in gray (solid = activation, dashed = inhibition)
    _draw_network_on_ax(axes[0], G_orig, pos_orig, orig_cmap,
                        activation_color, inhibition_color, node_size, font_size,
                        gray_edge_set=orig_edges,
                        missing_node_set=missing_nodes,
                        missing_edge_set=missing_edges)
    axes[0].set_title("Original Network", fontsize=12, fontweight='bold', pad=8)
    axes[0].axis('off')

    # Right panel: original edges gray, new KG edges colored and thicker
    _draw_network_on_ax(axes[1], G_ext, pos_shared, ext_cmap,
                        activation_color, inhibition_color, node_size, font_size,
                        new_node_set=new_nodes, new_edge_set=new_edges,
                        gray_edge_set=orig_edges, edge_scores=edge_scores)
    axes[1].set_title("Extended Network", fontsize=12, fontweight='bold', pad=8)
    axes[1].axis('off')

    # --- Per-panel legends ---
    # Left panel: node colours + edge styles
    left_handles = [Patch(facecolor=color_node, edgecolor='#333333', label='Original nodes')]
    for group_name, (_, gc) in node_groups.items():
        left_handles.append(Patch(facecolor=gc, edgecolor='#333333', label=group_name))
    left_handles += [
        Line2D([0], [0], color='#AAAAAA', lw=1.5, ls='-',  label='Activation'),
        Line2D([0], [0], color='#AAAAAA', lw=1.5, ls='--', label='Inhibition'),
    ]
    if missing_nodes:
        left_handles.append(Patch(facecolor=color_node + '99', edgecolor='#AAAAAA',
                                  linestyle='dashed', label='Missing nodes'))
    if missing_edges:
        left_handles.append(Line2D([0], [0], color='#CCCCCC', lw=1.5, ls=':', label='Missing edges'))
    axes[0].legend(handles=left_handles, loc='upper left', bbox_to_anchor=(0, 1),
                   fontsize=8, framealpha=0.85, borderaxespad=0.5)

    # Right panel: original + new node colours + edge styles distinguishing new vs existing
    right_handles = [
        Patch(facecolor=color_node,           edgecolor='#333333',             label='Original nodes'),
        Patch(facecolor='#FFFFFF', edgecolor='#FFC000', linewidth=2, label='New KG nodes'),
    ]
    for group_name, (_, gc) in node_groups.items():
        right_handles.append(Patch(facecolor=gc, edgecolor='#333333', label=group_name))
    right_handles += [
        Line2D([0], [0], color='#AAAAAA',        lw=1.5, ls='-',  label='Existing activation'),
        Line2D([0], [0], color='#AAAAAA',        lw=1.5, ls='--', label='Existing inhibition'),
        Line2D([0], [0], color=activation_color, lw=3.5, ls='-',  label='New activation'),
        Line2D([0], [0], color=inhibition_color, lw=3.5, ls='--', label='New inhibition'),
    ]
    axes[1].legend(handles=right_handles, loc='upper left', bbox_to_anchor=(0, 1),
                   fontsize=8, framealpha=0.85, borderaxespad=0.5)

    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=1.01)

    plt.tight_layout()

    if return_fig:
        return fig, axes
    plt.show()
    return None
