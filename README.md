# SIDH-search

## Installation

This library uses the reference imeplementation of SIKE submitted to NIST. Compile those libraries and copy them into a
folder with the name `lib` or in e.g. `/usr/local/lib`.

Both libraries require GMP.

## Graph automorphisms

For `graph_automorphisms.R` to be able to read the `.dot`-files produced from `explore_graph` it needs the arrow format
to be changed from `--` to the undirected `->`. This can be done by a simple application of `sed`.

``` {shell}
sed -i "" "s/--/->/" output/*
```
