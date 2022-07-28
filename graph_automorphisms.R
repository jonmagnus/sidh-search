
library("sna")
library("igraph")

# TODO: Handle big grpahs
args = commandArgs(trailingOnly=TRUE)
for (filename in args) {
    g.dot <- read.dot(filename)
    g.dot <- g.dot + t(g.dot)
    g.graph <- graph.adjacency(g.dot)
    print(filename)
    print(automorphisms(g.graph))
    auts = automorphism_group(g.graph)
    for (aut in auts) {
        cat("Permuted indices =",sum(aut != 1:length(aut)),"/",length(aut),"\n")
    }
}
