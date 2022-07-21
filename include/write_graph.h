#ifndef WRITE_GRAPH_H
#define WRITE_GRAPH_H

#include <stdio.h>
#include <fp2.h>

unsigned int hex_gradient(unsigned int start,
                          unsigned int end,
                          int max,
                          int val);

int pop_count(const int *mask, int n);

void filter_slope(const ff_Params *p,
                  const fp2 *slope,
                  fp2 *invariants,
                  int num_nodes,
                  int **edges,
                  int num_edges,
                  int *node_filter,
                  int *edge_filter);

void get_subgraph(const int *node_filter,
                  int num_nodes,
                  const int *edge_filter,
                  int  **edges,
                  int num_edges,
                  int *sub_nodes,
                  int *num_sub_nodes,
                  int **sub_edges,
                  int *num_sub_edges);

void write_graph(FILE *stream,
                 const int *depths,
                 char **keys,
                 const unsigned int *colors,
                 int num_nodes,
                 int **edges,
                 int num_edges);

#endif // WRITE_GRAPH_H
