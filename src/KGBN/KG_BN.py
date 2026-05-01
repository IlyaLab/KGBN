# building a boolean network

# given a KG subset, can we construct a boolean network?

# basic cache
GRAPH_TABLES = {}

def all_boolean_combos(n):
    "Generates all boolean strings of length n, as tuples of 1s and 0s."
    if n == 0:
        return []
    elif n == 1:
        return [(0,), (1,)]
    else:
        results = all_boolean_combos(n-1)
        new_results = []
        for r in results:
            new_results.append((0,) + r)
            new_results.append((1,) + r)
        return new_results


# try using signor
def load_signor_network(gene_list, input_format="symbol", joiner='&', kg_filename='SIGNOR_2025_08_14.tsv',
        only_proteins=True, score_cutoff=None, saveto=None):
    """
    Creates a boolean network from SigNOR using all of the provided genes. Tries to build a connected Steiner subgraph...

    Args:
        gene_list - list of gene symbols, gene ids, or uniprot ids.
        input_format - "symbol", "id", or "uniprot"
        joiner - "&", "|", "inhibitor_wins", or "majority", or "plurality" (difference between the last two: in plurality, a tie indicates 1, in majority, a tie indicates 0)
        kg_filename - "SIGNOR_2025_08_14.tsv" by default. "SIGNOR_formatted.tsv" can also be used (this is an older version of SigNOR).
        only_proteins - whether to only use protein nodes in SIGNOR for getting the subgraph (default: True)
        score_cutoff - minimum score threshold for edges to be included (default: None, includes all edges)
        saveto - if provided, save all edge attributes (source, target, direction, predicate,
                 Primary_Knowledge_Source, Knowledge_Source, publications, score, …) to this
                 CSV file path (e.g. 'KG_evidence.csv').  Default: None (no file written).
    """
    from . import graph_info, steiner_tree, gene_names
    if ' ' not in joiner and (joiner == '&' or joiner == '|'):
        joiner = ' ' + joiner + ' '
    input_format = input_format.lower()
    
    # Create a unique key for the graph table with score cutoff
    if score_cutoff is not None:
        graph_key = f"{kg_filename}_score_cutoff_{score_cutoff}"
    else:
        graph_key = kg_filename
    
    if graph_key in GRAPH_TABLES:
        graph_table = GRAPH_TABLES[graph_key]
    else:
        graph_table = graph_info.load_graph(kg_filename)
        
        # Apply score cutoff if specified
        if score_cutoff is not None:
            if 'score' in graph_table.columns:
                n_original = len(graph_table)
                graph_table = graph_table[graph_table['score'] >= score_cutoff].copy()
                print(f"Applied score cutoff {score_cutoff}, filtered to {len(graph_table)}/{n_original} edges")        
        GRAPH_TABLES[graph_key] = graph_table
    
    graph = graph_info.df_to_graph(graph_table, False)
    digraph = graph_info.df_to_graph(graph_table, True)
    if only_proteins:
        protein_names = [n['name'] for n in graph_info.nodes_in_category(graph, 'protein')]
        graph = graph.induced_subgraph(protein_names)
        digraph = digraph.induced_subgraph(protein_names)
    # get a graph subset
    # signor names are uniprot, of the format UNIPROT::[uniprot ID]
    # feature_name is the gene name
    uniprot_list = []
    if input_format == 'symbol':
        id_list = gene_names.get_ids(gene_list)
        print(f"number of genes found: {len(id_list)}")
        print(id_list)
        uniprot_list = gene_names.gene_ids_to_uniprot(id_list)
    elif input_format == 'id' or input_format == 'gene_id':
        uniprot_list = gene_names.gene_ids_to_uniprot(gene_list)
    else:
        uniprot_list = gene_list
    uniprot_list = ['UNIPROT::'+x for x in uniprot_list]
    # print(uniprot_list)
    # filter uniprot list based on existence in the graph
    new_uniprot_list = []
    for u in uniprot_list:
        try:
            graph.vs.find(u)
            new_uniprot_list.append(u)
        except:
            continue
    # # print edges in the graph where subject and object is in ['CDKN2A', 'MYC', 'TP53','CDK2']
    # print("\nEdge found among CDKN2A, MYC, TP53, CDK2:")
    # for e in graph.es:
    #     if graph.vs[e.source].attributes()['feature_name'] in ['CDKN2A', 'MYC', 'TP53','CDK2'] and graph.vs[e.target].attributes()['feature_name'] in ['CDKN2A', 'MYC', 'TP53','CDK2']:
    #         print(f"{graph.vs[e.source].attributes()['feature_name']} {e.attributes()['predicate']} {graph.vs[e.target].attributes()['feature_name']}")
    # print(new_uniprot_list)
    # Steiner subgraph - get a tree using a connected thing...
    tree = steiner_tree.steiner_tree(graph, new_uniprot_list)
    subgraph = digraph.induced_subgraph([n['name'] for n in tree.vs])
    # use the subgraph to build a connected thing
    # TODO: figure out the rule for combining inputs
    bn_lines = []
    # inhibitors
    # get all input edges
    all_relations = []
    all_edge_records = []  # for optional CSV export
    for n in subgraph.vs:
        inhibitors = []
        upregulators = []
        incoming_genes = set()
        gene_name = n.attributes()['feature_name']
        in_edges = subgraph.incident(n, mode='in')
        input_nodes = []
        edge_scores = []

        for e in in_edges:
            e = subgraph.es[e]
            in_node = subgraph.vs[e.source].attributes()['feature_name']
            if in_node in incoming_genes:
                continue
            incoming_genes.add(in_node)
            attrs = e.attributes()
            predicate = attrs.get('predicate', '')
            score = attrs.get('score', None)

            if 'down-regulates' in predicate:
                direction = 'inhibit'
                input_nodes.append(f'(! {in_node})')
                inhibitors.append(in_node)
                all_relations.append((in_node, gene_name, direction, score))
                if score is not None:
                    edge_scores.append(f"{in_node}_inhibit:{score}")
            elif 'up-regulates' in predicate:
                direction = 'activate'
                input_nodes.append(f'({in_node})')
                upregulators.append(in_node)
                all_relations.append((in_node, gene_name, direction, score))
                if score is not None:
                    edge_scores.append(f"{in_node}_activate:{score}")
            else:
                direction = 'unknown'

            # Collect full record for CSV
            record = {'source': in_node, 'target': gene_name, 'direction': direction}
            record.update(attrs)
            all_edge_records.append(record)
        
        if joiner.lower() in ['inhibitor_wins', 'inhibitorwins', 'inhibitor wins']:
            input_nodes_string = ''
            if len(inhibitors) > 1:
                inhibitor_string = ' & '.join(f'!{x}' for x in inhibitors)
                inhibitor_string =f'({inhibitor_string})'
            elif len(inhibitors) == 1:
                inhibitor_string = f'!{inhibitors[0]}'
            else:
                inhibitor_string = ''

            if len(upregulators) > 1:
                upregulator_string = ' | '.join(f'{x}' for x in upregulators)
                upregulator_string = f'({upregulator_string})'
            elif len(upregulators) == 1:
                upregulator_string = f'{upregulators[0]}'
            else:
                upregulator_string = ''

            if inhibitor_string and upregulator_string:
                input_nodes_string = f'{inhibitor_string} & {upregulator_string}'
            elif inhibitor_string:
                input_nodes_string = f'{inhibitor_string}'
            elif upregulator_string:
                input_nodes_string = f'{upregulator_string}'
            else:
                input_nodes_string = gene_name
        elif joiner.lower() == 'majority' or joiner.lower() == 'plurality':
            # a majority vote system for upregulators and downregulators
            # node is true if activators - repressors > 0
            input_nodes_string = ''
            up_down = upregulators + inhibitors
            all_outputs = []
            for permutation in all_boolean_combos(len(up_down)):
                upreg = sum(permutation[:len(upregulators)])
                downreg = sum(permutation[len(upregulators):])
                if (joiner == 'majority' and  upreg > downreg) or (joiner == 'plurality' and upreg >= downreg) :
                    gene_string = ' & '.join('!'+g if v==0 else g for v, g in zip(permutation, up_down))
                    all_outputs.append('('+gene_string+')')
            input_nodes_string = ' | '.join(all_outputs)
            if len(input_nodes_string) == 0:
                input_nodes_string = '0'
        else:
            input_nodes_string = joiner.join(input_nodes)

        if len(input_nodes) == 0:
            input_nodes_string = gene_name
        
        # Add score annotations as comments
        if edge_scores:
            score_comment = ' # Scores: ' + '; '.join(edge_scores)
            output_string = f'{gene_name} = {input_nodes_string}{score_comment}'
        else:
            output_string = f'{gene_name} = {input_nodes_string}'
        
        bn_lines.append(output_string)
    # order the equations by key alphabetically
    bn_lines.sort(key=lambda x: x.split('=')[0])

    if saveto and all_edge_records:
        import pandas as pd
        df = pd.DataFrame(all_edge_records)
        df.to_csv(saveto, index=False)
        print(f"Edge attributes saved to {saveto} ({len(df)} rows)")

    return '\n'.join(bn_lines), all_relations


def merge_PBN_string(original_string, KG_string, prob=0.5):
    """
    Merge the original model and the KG model to a PBN
    prob: probability of the equations from the original model
    """
    # Parse equations from both models
    original_equations = {}
    for line in original_string.strip().split('\n'):
        if '=' in line and not line.strip().startswith('#'):
            # Remove inline comments
            if '#' in line:
                line = line.split('#')[0].strip()
            target, rule = line.split('=', 1)
            original_equations[target.strip()] = rule.strip()
    
    kg_equations = {}
    for line in KG_string.strip().split('\n'):
        if '=' in line and not line.strip().startswith('#'):
            # Remove inline comments
            if '#' in line:
                line = line.split('#')[0].strip()
            target, rule = line.split('=', 1)
            kg_equations[target.strip()] = rule.strip()
    
    # Merge equations
    merged_equations = []
    all_targets = set(original_equations.keys()) | set(kg_equations.keys())
    
    for target in all_targets:
        if target in original_equations and target in kg_equations:
            if original_equations[target] == kg_equations[target]:
                # Both models have the same equation for this target
                merged_equations.append(f"{target} = {original_equations[target]}, 1")
            else:
                # Use both with specified probabilities
                merged_equations.append(f"{target} = {original_equations[target]}, {prob}")
                merged_equations.append(f"{target} = {kg_equations[target]}, {1-prob}")
        elif target in original_equations:
            # Only original model has this target
            merged_equations.append(f"{target} = {original_equations[target]}, 1")
            print(f"Only original model has this target: {target}")
        else:
            # Only KG model has this target
            print("Only KG model has this target:", target)
            merged_equations.append(f"{target} = {kg_equations[target]}, 1")
    
    # remove equation with prob = 0
    merged_equations = [eq for eq in merged_equations if eq.split(',')[1] != ' 0']
    merged_string = '\n'.join(merged_equations)
    return merged_string

def _get_edge_direction(predicate, qualifiers):
    """Infer activation/inhibition from a Translator edge's predicate and qualifiers.

    Returns 'inhibit', 'activate', or 'unknown'.
    """
    for q in (qualifiers or []):
        if q.get('qualifier_type_id') == 'biolink:object_direction_qualifier':
            val = q.get('qualifier_value', '')
            if val in ('downregulated', 'decreased'):
                return 'inhibit'
            if val in ('upregulated', 'increased'):
                return 'activate'
    pred = (predicate or '').lower()
    if 'negatively' in pred or 'negative' in pred:
        return 'inhibit'
    if 'positively' in pred or 'positive' in pred:
        return 'activate'
    return 'unknown'


def _resolve_curie_names(curies, curie_to_symbol, filtered_results, gene_names_mod, name_resolver):
    """Populate *curie_to_symbol* for every CURIE in *curies* that is not yet mapped."""
    # 1. Extract subject_name / object_name from edge attributes
    for v in filtered_results.values():
        for attr in v.get('attributes', []):
            tid = attr.get('attribute_type_id')
            if tid == 'subject_name' and v['subject'] not in curie_to_symbol:
                curie_to_symbol[v['subject']] = attr['value']
            elif tid == 'object_name' and v['object'] not in curie_to_symbol:
                curie_to_symbol[v['object']] = attr['value']

    # 2. NCBIGene CURIEs → gene symbols via local gene_names module
    for curie in curies:
        if curie not in curie_to_symbol and curie.startswith('NCBIGene:'):
            try:
                curie_to_symbol[curie] = gene_names_mod.get_symbol(int(curie.split(':')[1]))
            except (KeyError, ValueError):
                pass

    # 3. Remaining unknowns → Translator name_resolver
    unknown = [c for c in curies if c not in curie_to_symbol]
    if unknown:
        try:
            extra = name_resolver.batch_lookup(unknown)
            for c in unknown:
                info = extra.get(c)
                if info:
                    name = getattr(info, 'label', None) or getattr(info, 'name', None)
                    if name:
                        curie_to_symbol[c] = name
                        continue
                curie_to_symbol.setdefault(c, c)
        except Exception:
            for c in unknown:
                curie_to_symbol.setdefault(c, c)


def load_translator_network(gene_list, type=None, joiner='&', n_hop=1,
        sele_predicate=['biolink:regulates'], sele_tax_id=9606, sele_API=None, sele_KG=None):
    """
    Creates a boolean network from Translator using all of the provided genes.
    Tries to build a connected Steiner subgraph among input nodes (via n_hop paths).

    Args:
        gene_list - list of symbols for genes, proteins, chemicals, or drugs
        type - "gene", "protein", "chemical", "drug", "disease", "phenotype", or None (all possible types)
        joiner - "&", "|", "inhibitor_wins", or "majority", or "plurality"
                 (difference between the last two: in plurality, a tie indicates 1,
                 in majority, a tie indicates 0)
        n_hop - number of hops to build the subgraph, e.g., n_hop=1 means only
                direct connections between input, n_hop=2 means one hop away, etc.
                * currently only n_hop=1 is supported in TCT
        sele_predicate - list of predicates to query from Translator
        sele_tax_id - tax id to use when looking up genes/proteins, default is 9606 (human)
        sele_API - list of APIs to select from Translator
        sele_KG - list of KG to select from results (matches sources resource_id
                  or attribute provided_by / knowledge_source values)

    Returns:
        (bn_string, all_relations) where bn_string is a boolean network string
        and all_relations is a list of (source, target, direction, score) tuples.
    """
    from TCT import name_resolver, translator_metakg, translator_kpinfo, translator_query, TCT
    import igraph as ig
    from . import gene_names, steiner_tree

    if ' ' not in joiner and joiner in ('&', '|'):
        joiner = ' ' + joiner + ' '

    # 1. Resolve gene symbols → CURIEs
    print("Resolving gene names...")
    if sele_tax_id:
        sele_tax_id = "NCBITaxon:" + str(sele_tax_id)
        print(f"Selected tax id: {sele_tax_id}")
    else:
        print("No tax id selected, using all")
    input_info = name_resolver.batch_lookup(gene_list, only_taxa = sele_tax_id)
    input_curies = []
    curie_to_symbol = {}
    for g in gene_list:
        info = input_info.get(g)
        print(f"info for {g}: {info}")
        if info is not None and hasattr(info, 'curie') and info.curie:
            input_curies.append(info.curie)
            curie_to_symbol[info.curie] = g
    print(f"Resolved {len(input_curies)}/{len(gene_list)} genes to CURIEs")
    if len(input_curies) < 2:
        print("Need at least 2 resolved genes to build a network.")
        return '', []
    
    # determine the input node categories
    if type:
        category_map = {
            'gene': 'biolink:Gene',
            'protein': 'biolink:Protein',
            'chemical': 'biolink:ChemicalSubstance',
            'drug': 'biolink:Drug',
            'disease': 'biolink:DiseaseOrPhenotypicFeature',
            'phenotype': 'biolink:PhenotypicFeature',
        }
        input_node_category = [category_map.get(type)]
        if input_node_category is None:
            raise ValueError(f"Unknown type '{type}'. Must be one of: {list(category_map.keys())}")
    else:
        input_node_category = []
        for node in gene_list:
            if hasattr(input_info[node], 'types'):
                input_node_category.extend(input_info[node].types)
        input_node_category = list(set(input_node_category))
        print(f"Querying for input node categories: {input_node_category}")

    # 2. Load Translator resources & determine available predicates/APIs
    print("Loading Translator resources...")
    _, APInames = translator_kpinfo.get_translator_kp_info()
    metaKG = translator_metakg.get_KP_metadata(APInames)
    APInames, metaKG = translator_metakg.add_plover_API(APInames, metaKG)

    all_predicates = list(set(
        TCT.select_concept(sub_list=input_node_category, obj_list=input_node_category, metaKG=metaKG)))
    query_predicates = all_predicates
    if sele_predicate is not None:
        filtered_predicates = [p for p in all_predicates if p in sele_predicate]
        if filtered_predicates:
            query_predicates = filtered_predicates
        else:
            print(f"Warning: none of sele_predicate found. Using all {len(all_predicates)} available predicates.")
    
    all_APIs = TCT.select_API(sub_list=input_node_category, obj_list=input_node_category, metaKG=metaKG)
    query_APIs = all_APIs
    if sele_API is not None:
        filtered_APIs = [a for a in all_APIs if a in sele_API]
        if filtered_APIs:
            query_APIs = filtered_APIs
        else:
            print(f"Warning: none of sele_API found. Using all {len(all_APIs)} available APIs.")

    API_predicates = {}
    for api in set(metaKG['API']):
        API_predicates[api] = list(set(metaKG[metaKG['API'] == api]['Predicate']))

    # 3. Query Translator
    def _do_query(curies):
        qj = TCT.format_query_json(curies, [], input_node_category, input_node_category, query_predicates)
        return translator_query.parallel_api_query(
            query_json=qj, select_APIs=query_APIs,
            APInames=APInames, API_predicates=API_predicates,
            max_workers=len(query_APIs))

    print(f"Querying Translator with {len(input_curies)} input nodes...")
    all_results = _do_query(input_curies)
    print(f"Query returned {len(all_results)} results")

    # 4. Iterative expansion for n_hop > 2
    queried = set(input_curies)
    for hop_i in range(n_hop - 2):
        boundary = set()
        for v in all_results.values():
            if not isinstance(v, dict):
                continue
            s, o = v.get('subject'), v.get('object')
            if s in queried and o and o not in queried:
                boundary.add(o)
            if o in queried and s and s not in queried:
                boundary.add(s)
        if not boundary:
            print(f"No new boundary nodes at hop {hop_i + 3}.")
            break
        print(f"Hop {hop_i + 3}: querying {len(boundary)} boundary nodes...")
        new_results = _do_query(list(boundary))
        all_results.update(new_results)
        queried.update(boundary)
        print(f"Total results so far: {len(all_results)}")

    # 5. Filter results by KG
    sele_KG_set = set(sele_KG) if sele_KG else None

    filtered = {}
    for k, v in all_results.items():
        if not isinstance(v, dict):
            continue
        if sele_KG_set:
            match = False
            for src in v.get('sources', []):
                rid = src.get('resource_id', '')
                if rid in sele_KG_set or rid.split(':')[-1] in sele_KG_set:
                    match = True
                    break
            if not match:
                for attr in v.get('attributes', []):
                    if attr.get('attribute_type_id') in ('provided_by', 'knowledge_source'):
                        if attr.get('value') in sele_KG_set:
                            match = True
                            break
            if not match:
                continue
        filtered[k] = v

    # n_hop = 1: keep only edges between input nodes
    input_curie_set = set(input_curies)
    if n_hop == 1:
        filtered = {k: v for k, v in filtered.items()
                    if v.get('subject') in input_curie_set
                    and v.get('object') in input_curie_set}

    print(f"Filtered to {len(filtered)} edges")
    if not filtered:
        print("No edges found after filtering.")
        return '', []

    # 6. Resolve CURIEs → gene symbols
    all_curies = set()
    for v in filtered.values():
        all_curies.add(v['subject'])
        all_curies.add(v['object'])
    _resolve_curie_names(all_curies, curie_to_symbol, filtered, gene_names, name_resolver)

    # 7. Collect directed edges as (source_sym, target_sym, pred, quals)
    edges_data = []
    for k, v in filtered.items():
        s, o = v['subject'], v['object']
        if s == o:
            continue
        sn = curie_to_symbol.get(s, s)
        on = curie_to_symbol.get(o, o)
        edges_data.append((sn, on, v.get('predicate', ''), v.get('qualifiers', [])))
    if not edges_data:
        print("No valid directed edges found.")
        return '', []

    # 8. For n_hop > 1: build undirected graph, find Steiner tree
    if n_hop > 1:
        all_names = sorted({s for s, *_ in edges_data} | {o for _, o, *_ in edges_data})
        G = ig.Graph(directed=False)
        for name in all_names:
            G.add_vertex(name)
        seen_undirected = set()
        for s, o, *_ in edges_data:
            ek = tuple(sorted([s, o]))
            if ek not in seen_undirected:
                seen_undirected.add(ek)
                G.add_edge(s, o)

        input_names = []
        for curie in input_curies:
            sym = curie_to_symbol.get(curie)
            if sym:
                try:
                    G.vs.find(name=sym)
                    if sym not in input_names:
                        input_names.append(sym)
                except ValueError:
                    pass

        if len(input_names) >= 2:
            print(f"Finding Steiner tree for {len(input_names)} input nodes "
                  f"in a graph of {len(G.vs)} nodes...")
            tree = steiner_tree.steiner_tree(G, input_names)
            tree_nodes = {v['name'] for v in tree.vs}
            edges_data = [(s, o, p, q) for s, o, p, q in edges_data
                          if s in tree_nodes and o in tree_nodes]
            print(f"Steiner tree: {len(tree_nodes)} nodes, {len(edges_data)} directed edges")

    # 9. Group edges by target, determine direction, build BN equations
    # node_inputs: target_name → {source_name: [directions]}
    node_inputs = {}
    for sn, on, pred, quals in edges_data:
        direction = _get_edge_direction(pred, quals)
        if direction == 'unknown':
            print(f"Unknown direction for edge {sn} -> {on} with predicate {pred} and qualifiers {quals}")
            continue
        node_inputs.setdefault(on, {}).setdefault(sn, []).append(direction)

    all_node_names = set()
    for s, o, *_ in edges_data:
        all_node_names.add(s)
        all_node_names.add(o)

    bn_lines = []
    all_relations = []
    for gene_name in sorted(all_node_names):
        sources = node_inputs.get(gene_name, {})
        inhibitors = []
        upregulators = []
        input_nodes = []

        for src in sorted(sources):
            if src == gene_name:
                continue
            directions = sources[src]
            n_inh = sum(1 for d in directions if d == 'inhibit')
            n_act = len(directions) - n_inh
            if n_inh > n_act:
                input_nodes.append(f'(! {src})')
                inhibitors.append(src)
                all_relations.append((src, gene_name, 'inhibit', None))
            else:
                input_nodes.append(f'({src})')
                upregulators.append(src)
                all_relations.append((src, gene_name, 'activate', None))

        # --- build the equation string (same joiner logic as load_signor_network) ---
        if joiner.lower() in ['inhibitor_wins', 'inhibitorwins', 'inhibitor wins']:
            if len(inhibitors) > 1:
                inhibitor_string = ' & '.join(f'!{x}' for x in inhibitors)
                inhibitor_string = f'({inhibitor_string})'
            elif len(inhibitors) == 1:
                inhibitor_string = f'!{inhibitors[0]}'
            else:
                inhibitor_string = ''

            if len(upregulators) > 1:
                upregulator_string = ' | '.join(f'{x}' for x in upregulators)
                upregulator_string = f'({upregulator_string})'
            elif len(upregulators) == 1:
                upregulator_string = f'{upregulators[0]}'
            else:
                upregulator_string = ''

            if inhibitor_string and upregulator_string:
                input_nodes_string = f'{inhibitor_string} & {upregulator_string}'
            elif inhibitor_string:
                input_nodes_string = inhibitor_string
            elif upregulator_string:
                input_nodes_string = upregulator_string
            else:
                input_nodes_string = gene_name

        elif joiner.lower() in ['majority', 'plurality']:
            up_down = upregulators + inhibitors
            all_outputs = []
            for permutation in all_boolean_combos(len(up_down)):
                upreg = sum(permutation[:len(upregulators)])
                downreg = sum(permutation[len(upregulators):])
                if ((joiner.lower() == 'majority' and upreg > downreg)
                        or (joiner.lower() == 'plurality' and upreg >= downreg)):
                    gene_string = ' & '.join(
                        '!' + g if v == 0 else g
                        for v, g in zip(permutation, up_down))
                    all_outputs.append('(' + gene_string + ')')
            input_nodes_string = ' | '.join(all_outputs)
            if len(input_nodes_string) == 0:
                input_nodes_string = '0'
        else:
            input_nodes_string = joiner.join(input_nodes)

        if len(input_nodes) == 0:
            input_nodes_string = gene_name

        bn_lines.append(f'{gene_name} = {input_nodes_string}')

    bn_lines.sort(key=lambda x: x.split('=')[0])
    return '\n'.join(bn_lines), all_relations