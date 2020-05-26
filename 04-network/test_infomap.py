import networkx as nx
import infomap

def compute_infomap(G, name='net', path='results', markov_time=1):
    print('-- Compute Infomap --')
    myInfomap = infomap.Infomap("--two-level --markov-time 1 --silent --seed 123")
    # Map node-id to integer
    dict_node_to = {n: i for i, n in enumerate(G.nodes(), start=1)}
    dict_node_from = {i: n for n, i in dict_node_to.items()}
    #
    for i, d in G.nodes(data=True):
        it = dict_node_to[i]
        myInfomap.add_node(node_id=it, name=i)
        #myInfomap.set_name(node_id, name='')
        #myInfomap.add_state_node(state_id=1, node_id=it)
        #myInfomap.add_state_node(state_id=2, node_id=it)
        #myInfomap.add_state_node(state_id=2, node_id=it)

    for i, j, d in G.edges(data=True):
        # Map Translation
        it = dict_node_to[i]
        jt = dict_node_to[j]
        # from, to, weight
        myInfomap.add_link(source_id=it, target_id=jt, weight=d['weight'])

    # Run!
    myInfomap.run()
    # Debug print
    
    print("codelength:", myInfomap.codelength)
    print("num_leaf_modules:", myInfomap.num_leaf_modules)
    print("num_top_modules:", myInfomap.num_top_modules)
    print('tree')
    tree = myInfomap.tree
    for i in tree:
        print("child_index:",i.child_index, "depth:", i.depth, "module:", i.module_id, "path:", i.path)


    print('physical_tree')
    physical_tree = myInfomap.physical_tree
    for i in physical_tree:
        print("child_index:",i.child_index, "depth:", i.depth, "module:", i.module_id, "path:", i.path)

    print("get_nodes")
    nodes = myInfomap.get_nodes(depth_level=1, states=True)
    for i in nodes:
        print("child_index:",i.child_index, "depth:", i.depth, "module:", i.module_id, "path:", i.path)


    print("get_modules")
    modules = dict(myInfomap.get_modules(depth_level=0, states=True))
    print(modules)

    print("get_multilevel_modules")
    mmodules = dict(myInfomap.get_multilevel_modules(states=True))
    print(mmodules)
    
    print("Found {modules:d} modules with codelength: {codelength:.4f}".format(modules=myInfomap.num_top_modules, codelength=myInfomap.codelength))
    # Write Clu & Tree
    #myInfomap.writeMap("{path:s}/net_{name:s}.map".format(path=path, name=name))
    #myInfomap.writeClu("net_{name:s}.clu")
    #myInfomap.writeTree("{path:s}/net_{name:s}.tree"),
    # Dict of Results
    #dict_modules = myInfomap.get_modules(depth_level=1)
    #print(dict_modules)
    # translate nodes
    #dict_modules_translated = {dict_node_from[k]: v for k, v in dict_modules.items()}
    #
    return dict_modules_translated


if __name__ == '__main__':


    G = nx.Graph()

    G.add_nodes_from(['1', '2', '3'])
    G.add_nodes_from(['a', 'b', 'c'])

    #Intra
    G.add_edges_from([
        ('1', '2', {'weight':0.9}),
        ('1', '3', {'weight':0.9}),
        ('2', '3', {'weight':0.9})
    ])
    G.add_edges_from([
        ('a', 'b', {'weight':0.9}), 
        ('a', 'c', {'weight':0.9}), 
        ('b', 'c', {'weight':0.9})
    ])
    # Extra
    G.add_edges_from([
        ('a', '1', {'weight':0.5}), 
    ])
    

    print(G.nodes(data=True))
    print(G.edges(data=True))



    compute_infomap(G)
    