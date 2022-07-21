#include <stdio.h>
#include <stdlib.h>
#include <write_graph.h>
#include <fp2.h>

unsigned int hex_gradient(unsigned int start,
                          unsigned int end,
                          int max,
                          int val) {
    int rs, gs, bs;
    rs = (start >> 16) & 0xff;
    gs = (start >> 8) & 0xff;
    bs = start & 0xff;
    int re, ge, be;
    re = (end >> 16) & 0xff;
    ge = (end >> 8) & 0xff;
    be = end & 0xff;
    int r, g, b;
    r = val*(re - rs)/max + rs;
    g = val*(ge - gs)/max + gs;
    b = val*(be - bs)/max + bs;
    unsigned int h = ((r << 16) | (g << 8) | b) & 0xffffff;
    return h;
}

int pop_count(const int *mask, int n) {
    int pop = 0;
    for (int i = 0; i < n; i++) {
        if (mask[i]) pop++;
    }
    return pop;
}

/** 
 * Write the subgraph corresponding to curves with j-invariant fp-proporsional
 * to slope.
 */
void filter_slope(const ff_Params *p,
                  const fp2 *slope,
                  fp2 *invariants,
                  int num_nodes,
                  int **edges,
                  int num_edges,
                  int *node_filter,
                  int *edge_filter) {
    mp t1, t2;
    fp_Init(p, t1);
    fp_Init(p, t2);
    
    for (int i = 0; i < num_nodes; i++) {
        fp_Multiply(p, slope->x0, invariants[i].x1, t1);
        fp_Multiply(p, slope->x1, invariants[i].x0, t2);
        node_filter[i] = fp_IsEqual(p, t1, t2);
    }

    for (int i = 0; i < num_edges; i++) {
        int u, v;
        u = edges[i][0];
        v = edges[i][1];
        edge_filter[i] = node_filter[u] && node_filter[v];
    }

    fp_Clear(p, t1);
    fp_Clear(p, t2);
}

void get_subgraph(const int *node_filter,
                  int num_nodes,
                  const int *edge_filter,
                  int  **edges,
                  int num_edges,
                  int *sub_nodes,
                  int *num_sub_nodes,
                  int **sub_edges,
                  int *num_sub_edges) {
    int *reverse_indices = malloc(num_nodes*sizeof(int));
    *num_sub_nodes = 0;
    *num_sub_edges = 0;
    for (int i = 0; i < num_nodes; i++) {
        if (node_filter[i]) {
            reverse_indices[i] = *num_sub_nodes;
            sub_nodes[*num_sub_nodes] = i;
            (*num_sub_nodes)++;
        }
    }
    for (int i = 0; i < num_edges; i++) {
        if (edge_filter[i]) {
            int u, v, u_, v_;
            u = edges[i][0];
            v = edges[i][1];
            u_ = reverse_indices[u];
            v_ = reverse_indices[v];
            sub_edges[*num_sub_edges][0] = u_;
            sub_edges[*num_sub_edges][1] = v_;
            (*num_sub_edges)++;
        }
    }
    free(reverse_indices);
}

void write_graph(FILE *stream,
                 const int *depths,
                 char **keys,
                 const unsigned int *colors,
                 int num_nodes,
                 int **edges,
                 int num_edges) {
    int max_depth = 0;
    for (int i = 0; i < num_nodes; i++) {
        max_depth = depths[i] > max_depth ? depths[i] : max_depth;
    }
    fprintf(stream, "graph {\n");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(stream,
                "n%d [label=\"%s\\n%d\", color=\"#%06x\", style=filled]\n",
                i,
                keys[i],
                depths[i],
                colors[i]);
    }
    for (int i = 0; i < num_edges; i++) {
        fprintf(stream, "n%d--n%d\n", edges[i][0], edges[i][1]);
    }
    fprintf(stream, "}\n");
}
