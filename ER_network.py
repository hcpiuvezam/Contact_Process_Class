import networkx as nx
import numpy as np

N = 10 ** 6
p = 20 / N
er = nx.erdos_renyi_graph(N, p)

N_neigh = []
for key, value in er.adj.items():
    N_neigh.append(len(value))

print(np.mean(N_neigh))

indices_mat = []
for key, value in er.adj.items():
    indices_vec = []
    for i in value.keys():
        indices_vec.append(i)
    indices_mat.append(indices_vec)

print(len(indices_mat))

f = open("Networks/adj_ind_ER_10_6__20.txt", "w")
f.write('\n'.join([''.join([' {}'.format(item) for item in row]) for row in indices_mat]))
f.close()