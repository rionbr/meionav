

import networkx as nx

G = nx.Graph()

G.add_nodes_from(['1', '2', '3'], **{'mean-fert-rate': 0.5})
G.add_nodes_from(['a', 'b', 'c'])

G.add_edges_from([('1', '2'), ('1', '3'), ('2', '3')], **{'type': 'intra'})
G.add_edges_from([('a', 'b'), ('a', 'c'), ('b', 'c')], **{'type': 'intra'})
G.add_edges_from([('1', 'a'), ('2', 'b'), ('3', 'c')], **{'type': 'cross'})

print(G.nodes(data=True))
print(G.edges(data=True))
print(G.nodes['a'])


def transpose_variable_across_layers(G, variable):
    print('Transposing DM fertility results to HS & MM layers')
    dict_i_values = {i: d['mean-fert-rate'] for i, d in G.nodes(data=True) if d.get(variable, None) is not None}
    dict_j_values = {}
    for i, v in dict_i_values.items():
        cross_edges = [j for _i, j, d in G.edges(i, data=True) if d.get('type', None) == 'cross']
        for j in cross_edges:
            if j in dict_j_values:
                raise TypeError("Error due to multiple cross layer attribution of property. Multiple cross edges related to the same node")
            else:
                dict_j_values[j] = v
    nx.set_node_attributes(G, values=dict_j_values, name=variable)
    return G

G = transpose_variable_across_layers(G, 'mean-fert-rate')
print(G.nodes['a'])