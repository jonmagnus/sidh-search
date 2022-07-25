#include <gmp.h>
#include <search.h>
#include <sike_params_small.h>
#include <montgomery.h>
#include <isogeny.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mont_utils.h>
#include <printing.h>
#include <get_initial_curve.h>
#include <write_graph.h>

#define NUM_NODES 0xffffffllu
#define H_KEY_SIZE 40
#define H_KEY_OFFSET 20
#define EDGE_LABEL_SIZE 40
#define EDGE_STR_SIZE (2*H_KEY_SIZE + EDGE_LABEL_SIZE + 10)

// TODO: Speed up by avoiding sprintf.
// TODO: Implement heuristics for graph comparison.
// TODO: Implement search in 3-graph.
// TODO: Traverse graph by explicitly finding roots.
// TODO: Read $IKEp182 for inspiration on fast graph exploration.
// TODO: Embed on a periodic space (e.g. torus or finite field euclidean grid)
// in hope of graph folding/filtration properties.
// TODO: Find what choices of origin on a genus 1 curve give supersingular curves
// -- what is the supersingular locus of the moduli space of elliptic curves.

void fp2_get_key(const fp2 *raw_key, char *hkey) {
    sprintf(hkey,
            "%s+i*%s",
            mpz_get_str(NULL, 0, raw_key->x0),
            mpz_get_str(NULL, 0, raw_key->x1));
}

/**
 * Perform a BF-search to map the whole isogeny graph.
 * The isogeny graph over Fp2 should be connected,
 * so it will uncover the whole graph.
 */
int isogeny_bfs(const sike_params_t *params,
                mont_curve_int_t **q,
                int *end,
                int **edges,
                unsigned long long *num_edges,
                const int is_alice) {
    int (*reduce_to_single_torsion)(const mont_curve_int_t*,
                                    const mont_pt_t*,
                                    mont_pt_t*) = (
                                        is_alice
                                        ? &reduce_to_2_torsion
                                        : &reduce_to_3_torsion);
    ff_Params *p = params->EA.ffData;
    int rc = 0, front = 0;
    fp2 raw_key = { 0 };
    fp2_Init(p, &raw_key);
    char key[H_KEY_SIZE];
    ENTRY item = { .key = key };

    mont_curve_int_t *curve = NULL, *target_curve;

    mont_pt_t P_ = { 0 }, Q_ = { 0 }, PQ_ = { 0 }, PPQ_ = { 0 }, T={ 0 };
    mont_pt_t* iso_points[4] = { &P_, &Q_, &PQ_, &PPQ_ };
    for (int i = 0; i < 4 - is_alice; i++) mont_pt_init(p, iso_points[i]);
    mont_pt_init(p, &T);

    while(front < *end) {
        if (front % 100 == 0) {
            fprintf(stderr, "%08d/%08d wavefront %08d\r", front, *end, *end - front);
            fflush(stderr);
        }
        curve = q[front];
        j_inv(p, curve, &raw_key);
        item.key = malloc(H_KEY_SIZE);
        fp2_get_key(&raw_key, item.key);
        
        // Populate iso_points with {2,3}-torsion elements.
        ENTRY *found_item = hsearch(item , FIND);
        if (found_item == NULL) {
            fprintf(stderr, "Couldn't find item corresponding to front of queue\n");
            rc = 1;
            goto cleanup;
        }
        unsigned long long *data = found_item->data;
        free(item.key);
        int ordP, ordQ;
        ordP = (*reduce_to_single_torsion)(curve, &curve->P, &P_);
        ordQ = (*reduce_to_single_torsion)(curve, &curve->Q, &Q_);

        if (ordP < 2) {
            find_basis(curve, params->eA, params->eB, is_alice, &curve->P);
            ordP = (*reduce_to_single_torsion)(curve, &curve->P, &P_);
        }
        if (ordQ < 2) {
            find_basis(curve, params->eA, params->eB, is_alice, &curve->Q);
            ordQ = (*reduce_to_single_torsion)(curve, &curve->Q, &Q_);
        }

        if (is_alice) {
            // TODO: This need not cover all isogenies if P, Q generate same isogeny.
            while (fp2_IsEqual(p, &P_.x, &Q_.x)) {
                if (ordP > ordQ) {
                    xDBLe(curve, &curve->P, ordP - ordQ, &T);   // T and Q have the same order.
                    if (fp2_IsEqual(p, &curve->Q.x, &T.x)) {
                        find_basis(curve, params->eA, params->eB, is_alice, &curve->Q);
                        ordQ = (*reduce_to_single_torsion)(curve, &curve->Q, &Q_);
                    } else {
                        xADD(curve, &T, &curve->Q, &curve->Q);
                        ordQ = (*reduce_to_single_torsion)(curve, &curve->Q, &Q_);
                    }
                } else {
                    xDBLe(curve, &curve->Q, ordQ - ordP, &T);   // T and P have the same order.
                    if (fp2_IsEqual(p, &curve->P.x, &T.x)) {
                        find_basis(curve, params->eA, params->eB, is_alice, &curve->P);
                        ordP = (*reduce_to_single_torsion)(curve, &curve->P, &P_);
                    } else {
                        xADD(curve, &T, &curve->P, &curve->P);
                        ordP = (*reduce_to_single_torsion)(curve, &curve->P, &P_);
                    }
                }
            }

            if (fp2_IsConst(p, &P_.x, 0, 0)) {
                fp2_Invert(p, &Q_.x, &PQ_.x);
                //mont_pt_copy(p, &T, &curve->P);
            } else if (fp2_IsConst(p, &Q_.x, 0, 0)) {
                fp2_Invert(p, &P_.x, &PQ_.x);
                //mont_pt_copy(p, &T, &curve->Q);
            } else {
                fp2_Set(p, &PQ_.x, 0, 0);
            }
        } else {
            // Make sure P and Q do not give the same isogeny.
            while (fp2_IsEqual(p, &P_.x, &Q_.x)) {
                if (fp2_IsEqual(p, &P_.y, &Q_.y)) {
                    fp2_Negative(p, &curve->P.y, &curve->P.y);
                }
                if (ordP > ordQ) {
                    xTPLe(curve, &curve->P, ordP - ordQ, &T);   // T and Q have the same order.
                    if (fp2_IsEqual(p, &curve->Q.x, &T.x)) {
                        find_basis(curve, params->eA, params->eB, is_alice, &curve->Q);
                        ordQ = (*reduce_to_single_torsion)(curve, &curve->Q, &Q_);
                    } else {
                        xADD(curve, &T, &curve->Q, &curve->Q);
                        ordQ = (*reduce_to_single_torsion)(curve, &curve->Q, &Q_);
                    }
                } else {
                    xTPLe(curve, &curve->Q, ordQ - ordP, &T);   // T and P have the same order.
                    if (fp2_IsEqual(p, &curve->P.x, &T.x)) {
                        find_basis(curve, params->eA, params->eB, is_alice, &curve->P);
                        ordP = (*reduce_to_single_torsion)(curve, &curve->P, &P_);
                    } else {
                        xADD(curve, &T, &curve->P, &curve->P);
                        ordP = (*reduce_to_single_torsion)(curve, &curve->P, &P_);
                    }
                }
            }
            xADD(curve, &P_, &Q_, &PQ_);
            xADD(curve, &P_, &PQ_, &PPQ_);
        }
        for (int i = 0; i < 4 - is_alice; i++) {
            // Check isogeny for each principal torsion element.
            const mont_pt_t *T = iso_points[i];
            if (fp2_IsConst(p, &T->x, 0, 0)) {
                // Skip the case when evaluating 2-isogeny of (0,0).
                continue;
            }
            target_curve = malloc(sizeof(mont_curve_int_t));
            mont_curve_init(p, target_curve);
            if (is_alice) {
                curve_2_iso(p, T, curve, target_curve);
                eval_2_iso(p, T, &curve->P, &target_curve->P);
                eval_2_iso(p, T, &curve->Q, &target_curve->Q);
            } else {
                curve_3_iso(p, T, curve, target_curve);
                eval_3_iso(p, T, &curve->P, &target_curve->P);
                eval_3_iso(p, T, &curve->Q, &target_curve->Q);
            }

            j_inv(p, target_curve, &raw_key);
            item.key = malloc(H_KEY_SIZE);
            fp2_get_key(&raw_key, item.key);
            unsigned long long *target_data;
            ENTRY *target_item = hsearch(item, FIND);
            if (target_item == NULL) {
                // Add curve to table and queue.
                target_data = malloc(2*sizeof(unsigned long long));
                target_data[0] = data[0] + 1;   // Search depth
                target_data[1] = *end;
                item.data = target_data;

                if(hsearch(item, ENTER) == NULL) {
                    fprintf(stderr, "Failed to set table value");
                    rc = 1;
                    goto cleanup;
                }
                edges[*num_edges][0] = ((unsigned long long*)found_item->data)[1];
                edges[*num_edges][1] = *end;
                (*num_edges)++;
                q[(*end)++] = target_curve;
            } else {
                target_data = target_item->data;
                if (target_data[0] > data[0]) {
                    edges[*num_edges][0] = ((unsigned long long*)found_item->data)[1];
                    edges[*num_edges][1] = target_data[1];
                    (*num_edges)++;
                }
                free(item.key);
                mont_curve_clear(p, target_curve);
            }
        }
        // Clean up so only wavefront is in memory.
        //mont_curve_clear(p, q[front]);
        front++;
    }
    
cleanup:
    //for (int i = front; i < *end; i++) mont_curve_clear(p, q[i]);
    for (int i = 0; i < 3; i++) mont_pt_clear(p, iso_points[i]);
    mont_pt_clear(p, &T);
    fp2_Clear(p, &raw_key);
    fprintf(stderr, "Finished bfs cleanup\n");
    
    return rc;
}

int main(int argc, char **argv) {
    int rc = 0, is_alice = 0;
    if (argc >= 3) {
        // is_alice is passed.
        is_alice = strtol(argv[2], NULL, 0);
    }
    unsigned long long num_edges = 0;
    char hkey[H_KEY_SIZE];
    int *edges_ = malloc(2*NUM_NODES*sizeof(int));
    int **edges = malloc(NUM_NODES*sizeof(int*));
    if (edges_ == NULL || edges == NULL) {
        fprintf(stderr, "Failed to allocate edges\n");
        rc = 1;
        goto cleanup;
    }
    for (int i = 0; i < NUM_NODES; i++) {
        edges[i] = edges_ + 2*i;
    }
    fprintf(stderr, "Finished allocating edges");
    mont_curve_int_t **q = NULL;
    sike_params_raw_t raw_params = { 0 };
    sike_params_t params = { 0 };
    get_initial_curve(4, 3, &raw_params);
    //get_initial_curve(7, 9, &raw_params);
    mount_generic_bases(&raw_params);
    //sike_setup_params(&SIKEp33, &params);
    sike_setup_params(&raw_params, &params);

    rc = !hcreate(NUM_NODES);
    if (rc) {
        fprintf(stderr, "Failed to initialize table\n");
        goto cleanup;
    }

    q = malloc(NUM_NODES*sizeof(mont_curve_int_t*));
    memset(q, 0, NUM_NODES*sizeof(mont_curve_int_t*));
    
    mont_curve_int_t *initial_curve = (is_alice ? &params.EA : &params.EB);
    ff_Params *p = initial_curve->ffData;
    q[0] = initial_curve;

    fp2 raw_key = { 0 };
    fp2_Init(p, &raw_key);
    j_inv(p, initial_curve, &raw_key);
    fp2_get_key(&raw_key, hkey);
    unsigned long long initial_data[2] = { 0, 0 };
    ENTRY item = {
        .key = hkey,
        .data = initial_data,
    };
    hsearch(item, ENTER);

    int end = 1;
    isogeny_bfs(&params, q, &end, edges, &num_edges, is_alice);
    

    /////////////////
    // Parse result
    /////////////////

    fp2 *invariants = malloc(end*sizeof(fp2));
    int *depths = malloc(end*sizeof(int));
    char *keys_ = malloc(end*H_KEY_SIZE);
    char **keys = malloc(end*sizeof(char*));
    for (int i = 0; i < end; i++) {
        keys[i] = keys_ + i*H_KEY_SIZE;
    }
    for (int i = 0; i < end; i++) {
        const mont_curve_int_t *curve = q[i];
        j_inv(p, curve, &raw_key);
        fp2_Init(p, &invariants[i]);
        fp2_Copy(p, &raw_key, &invariants[i]);
        fp2_get_key(&raw_key, keys[i]);
        item.key = keys[i];
        ENTRY *found_item = hsearch(item, FIND);
        unsigned long long *data = found_item->data;
        depths[i] = (int)data[0];
    }
    
    ////////////////
    // Find subgraph
    ////////////////

    // Find optimal slope by pop-count.

    fp2 slope = { 0 };
    fp2_Init_set(p, &slope, 1, 0);
    int *node_filter, *edge_filter;
    node_filter = malloc(end*sizeof(int));
    edge_filter = malloc(num_edges*sizeof(int));
#if 0 // Turns out there are no other optimal linear embeddings.
    int best_pop = 1, best_pop_idx = 0;
    for (int i = 0; i < end; i++) {
        j_inv(p, q[i], &slope);
        if (fp_IsConstant(p, slope.x1, 0)) {
            // Ignore rational slopes.
            continue;
        }
        filter_slope(p, &slope, invariants, end, edges, num_edges, node_filter, edge_filter);
        int pop = pop_count(node_filter, end);
        if (pop > best_pop) {
            best_pop = pop;
            best_pop_idx = i;
        }
    }
    j_inv(p, q[best_pop_idx], &slope);
#endif
    filter_slope(p, &slope, invariants, end, edges, num_edges, node_filter, edge_filter);
    fp2_Clear(p, &slope);
    int num_sub_nodes, num_sub_edges;
    num_sub_nodes = pop_count(node_filter, end);
    num_sub_edges = pop_count(edge_filter, num_edges);
    int *sub_nodes, *sub_edges_, **sub_edges;
    sub_nodes = malloc(num_sub_nodes*sizeof(int));
    sub_edges_ = malloc(2*num_sub_edges*sizeof(int));
    sub_edges = malloc(num_sub_edges*sizeof(int*));
    for (int i = 0; i < num_sub_edges; i++) {
        sub_edges[i] = sub_edges_ + 2*i;
    }
    get_subgraph(node_filter,
                 end,
                 edge_filter,
                 edges,
                 num_edges,
                 sub_nodes,
                 &num_sub_nodes,
                 sub_edges,
                 &num_sub_edges);
        
    for (int i = 0; i < end; i++) {
        fp2_Clear(p, &invariants[i]);
    }
    free(invariants);

    //////////////
    // Write graph
    //////////////
    
    int *sub_depths = malloc(num_sub_nodes*sizeof(int));
    char **sub_keys = malloc(num_sub_nodes*sizeof(char*));
    for (int i = 0; i < num_sub_nodes; i++) {
        int u = sub_nodes[i];
        sub_depths[i] = depths[u];
        sub_keys[i] = keys[u];
    }

    FILE *logstream = stdout;
    if (argc >= 2) {
        logstream = fopen(argv[1], "w");
        if (logstream == NULL) {
            fprintf(stderr, "Failed to open file '%s'\n", argv[1]);
            rc = 1;
            goto cleanup;
        }
    }
    int max_depth = 0;
    for (int i = 0; i < end; i++) {
        max_depth = depths[i] > max_depth ? depths[i] : max_depth;
    }
    unsigned int *colors;
#if 0
    colors = malloc(end*sizeof(unsigned int));
    for (int i = 0; i < end; i++) {
        if (node_filter[i]) {
            colors[i] = 0x00ff00;
        } else {
            colors[i] = hex_gradient(0xff0000, 0x0000ff, max_depth, depths[i]);
        }
    }
    write_graph(logstream, depths, keys, colors, end, edges, num_edges);
#else
    colors = malloc(num_sub_nodes*sizeof(unsigned int));
    for (int i = 0; i < num_sub_nodes; i++) {
        colors[i] = hex_gradient(0xff0000, 0x0000ff, max_depth, sub_depths[i]);
    }
    write_graph(logstream,
                sub_depths,
                sub_keys,
                colors,
                num_sub_nodes,
                sub_edges,
                num_sub_edges);
#endif
    if (logstream != stdout) {
        fclose(logstream);
    }
    free(depths);
    free(keys);
    free(sub_keys);
    free(sub_depths);
    free(keys_);
    free(sub_nodes);
    free(sub_edges);
    free(sub_edges_);
    free(colors);
    free(node_filter);
    free(edge_filter);
    
cleanup:
    fp2_Clear(p, &raw_key);
    sike_teardown_params(&params);
    //hdestroy();
    free(q);
    free(edges);
    free(edges_);
    return rc;
}
