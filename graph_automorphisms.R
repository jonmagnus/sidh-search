
library("sna")
library("igraph")

args = commandArgs(trailingOnly=TRUE)
for (filename in args) {
    g.dot <- read.dot(filename)
    g.dot <- g.dot + t(g.dot)
    g.graph <- graph.adjacency(g.dot)
    print(automorphisms(g.graph))
    print(automorphism_group(g.graph))
}
