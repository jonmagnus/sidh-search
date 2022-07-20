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
            num_nodes += 1

        if b in node_map:
            v = node_map[b]
        else:
            v = node_map[b] = num_nodes
            labels.append(b)
            num_nodes += 1

        if u != v:
            adj_list.append([u, v])


adj_list = np.asarray(adj_list)
adj_list = np.sort(adj_list, axis=1)
adj_list = np.unique(adj_list, axis=0)
print(adj_list.shape)

adj_matrix = scipy.sparse.coo_matrix(
    (-np.ones_like(adj_list[:, 0]), (adj_list[:, 0], adj_list[:, 1])),
    shape=[num_nodes, num_nodes],
    dtype=float,
).tocsr()
adj_matrix = adj_matrix + adj_matrix.transpose()
adj_matrix.setdiag(-np.squeeze(np.asarray(adj_matrix.sum(0))))

dim = 3
w, v = scipy.sparse.linalg.eigsh(adj_matrix, dim + 1, which="LM", sigma=0)
print(w)
line_segments = v[:, 1:][adj_list]
fig = plt.figure()
if dim == 2:
    lc = mc.LineCollection(line_segments)
    ax = fig.add_subplot(111)
elif dim == 3:
    lc = Line3DCollection(line_segments)
    ax = fig.add_subplot(111, projection="3d")
else:
    raise f"Invalid dimension '{dim}'"

ax.plot(*(v[:, i] for i in range(1, dim + 1)), ".")
ax.add_collection(lc)
plt.show()
