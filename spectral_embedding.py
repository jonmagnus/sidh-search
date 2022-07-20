import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import sys

adj_list = []
node_map = {}
labels = []
is_rational = []
num_nodes = 0
with open(sys.argv[1], "r") as infile:
    infile.readline()
    for line in infile.readlines()[:-1]:
        [a, b] = line.strip().split("--")
        if a in node_map:
            u = node_map[a]
        else:
            u = node_map[a] = num_nodes
            labels.append(a)
            is_rational.append(b.split("i*")[1] == "0")
            num_nodes += 1

        if b in node_map:
            v = node_map[b]
        else:
            v = node_map[b] = num_nodes
            labels.append(b)
            is_rational.append(b.split("i*")[1] == '0"')
            num_nodes += 1

        if u != v:
            adj_list.append([u, v])


adj_list = np.asarray(adj_list)
adj_list = np.sort(adj_list, axis=1)
adj_list = np.unique(adj_list, axis=0)

adj_matrix = scipy.sparse.coo_matrix(
    (-np.ones_like(adj_list[:, 0]), (adj_list[:, 0], adj_list[:, 1])),
    shape=[num_nodes, num_nodes],
    dtype=float,
).tocsr()
adj_matrix = adj_matrix + adj_matrix.transpose()
adj_matrix.setdiag(-np.squeeze(np.asarray(adj_matrix.sum(0))))

dim = 3
w, v = scipy.sparse.linalg.eigsh(adj_matrix, dim + 1, which="LM", sigma=0)
is_rational = np.asarray(is_rational)
filter_adj_list = adj_list[is_rational[adj_list].all(1)]
# line_segments = v[:, 1:][adj_list]
line_segments = v[:, 1:][filter_adj_list]
fig = plt.figure()
if dim == 2:
    lc = mc.LineCollection(line_segments)
    ax = fig.add_subplot(111)
elif dim == 3:
    lc = Line3DCollection(line_segments)
    ax = fig.add_subplot(111, projection="3d")
else:
    raise f"Invalid dimension '{dim}'"

ax.scatter(*(v[is_rational, i] for i in range(1, dim + 1)), c="r")
ax.scatter(*(v[np.logical_not(is_rational), i] for i in range(1, dim + 1)), c="b")
ax.add_collection(lc)
plt.show()
