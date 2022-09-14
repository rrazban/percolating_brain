"""
Visualize random graph (nodes and edges)
at an average degree of 1. 

"""


import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

layout = nx.spring_layout


#n = 64  #corresponds to Harvard-Oxford atlas
#p_giant = 1.0 / (n )    #corresponds to k=1

n = 100
p_giant = 50.0 / (n )    #corresponds to k=1


for p in [p_giant]:
    plt.title(f"Random Graph Simulation", fontsize=16)

    G = nx.binomial_graph(n, p)
    if True: 
            print('average degree: {0:.2f}'.format(np.mean([deg for (node, deg) in G.degree()])))
            print('average clustering coefficient: {0:.2f}'.format(nx.average_clustering(G)))
            print('average shortest path length: {0}'.format(nx.average_shortest_path_length(G)))


    pos = layout(G)
    G.remove_nodes_from(list(nx.isolates(G)))	#remove orphans from visualization
    nx.draw(G, pos, with_labels=False, node_size=20)

    # identify largest connected component
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G0 = G.subgraph(Gcc[0])
    nx.draw_networkx_edges(G0, pos, edge_color="r", width=6.0)

    # show other connected components
    for Gi in Gcc[1:]:
        if len(Gi) > 1:
            nx.draw_networkx_edges(
                G.subgraph(Gi), pos, edge_color="r", alpha=0.3, width=5.0,
            )

plt.tight_layout()
plt.show()
