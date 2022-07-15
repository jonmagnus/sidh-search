#include <gmp.h>
#include <search.h>
#include <sike_params.h>
#include <sike_params_small.h>
#include <montgomery.h>
#include <isogeny.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mont_utils.h>
#include <printing.h>

#define NUM_NODES 0xfffff
#define H_KEY_SIZE 30
#define H_KEY_OFFSET 15

void fp2_get_key(const fp2 *raw_key, char *hkey) {
    sprintf(hkey,
            "%10s + i*%10s",
            mpz_get_str(NULL, 0, raw_key->x0),
            mpz_get_str(NULL, 0, raw_key->x1));
    //mpz_get_str(hkey, 16, raw_key->x0);
    //mpz_get_str(hkey + H_KEY_OFFSET, 16, raw_key->x1);
}

int isogeny_bfs(ff_Params *p,
                mont_curve_int_t **q,
                int end,
                char **edges,
                unsigned long long *num_edges) {
    int rc = 0, front = 0;
    /////////////////////////////////////////
    // Perform bfs and store edges as strings
    /////////////////////////////////////////
    fp2 raw_key = { 0 };
    fp2_Init(p, &raw_key);
    char key[H_KEY_SIZE];
    ENTRY item = { .key = key };

    mont_curve_int_t *curve = NULL, *target_curve;

    mont_pt_t P_ = { 0 }, Q_ = { 0 }, PQ_ = { 0 };
    mont_pt_t* iso2_points[3] = { &P_, &Q_, &PQ_ };
    for (int i = 0; i < 3; i++) mont_pt_init(p, iso2_points[i]);

    while (front < end && end < 10) {
        fprintf(stderr, "front=%5d, end=%5d\n", front, end);
        curve = q[front];
        printf("\n################\nFront curve\n################\n");
        mont_curve_printf(curve);
        printf("P: ");
        mont_pt_printf(&curve->P);
        printf("Q: ");
        mont_pt_printf(&curve->Q);
        j_inv(p, curve, &raw_key);
        fp2_get_key(&raw_key, item.key);
        fprintf(stderr, "key: %s\n", item.key);
        
        // Populate iso2_points with 2-torsion elements.
        ENTRY *found_item = hsearch(item , FIND);
        if (found_item == NULL) {
            fprintf(stderr, "Couldn't find item corresponding to front of queue\n");
            rc = 1;
            goto end;
        }
        int ordP, ordQ;
        ordP = reduce_to_2_torsion(curve, &curve->P, &P_);
        ordQ = reduce_to_2_torsion(curve, &curve->Q, &Q_);
        fprintf(stderr, "ordP=%d, ordQ=%d\n", ordP, ordQ);
        xADD(curve, &P_, &Q_, &PQ_);

        for (int i = 0; i < 3; i++) {
            printf("\nis2_points %d: ", i);
            mont_pt_printf(iso2_points[i]);
            if (!mont_is_principal_2_torsion(curve, iso2_points[i])){
                fprintf(stderr, "Point is not principal 2-torsion\n");
                rc = 1;
                goto end;
            }
            if (!mont_is_on_curve(curve, iso2_points[i])) {
                printf("iso2 point not on curve\n");
            }
            // Check the three possible 2-torsion elements.
            const mont_pt_t *T = iso2_points[i];
            target_curve = malloc(sizeof(mont_curve_int_t));
            mont_curve_init(p, target_curve);
            curve_2_iso(p, T, curve, target_curve);
            eval_2_iso(p, T, &curve->P, &target_curve->P);
            eval_2_iso(p, T, &curve->Q, &target_curve->Q);
            if (!mont_is_on_curve(target_curve, &target_curve->P)) {
                printf("Target P not on target curve: ");
                mont_pt_printf(&target_curve->P);
            }
            if (!mont_is_on_curve(target_curve, &target_curve->Q)) {
                printf("Target Q not on target curve: ");
                mont_pt_printf(&target_curve->Q);
            }

            int t_ordP, t_ordQ;
            t_ordP = reduce_to_2_torsion(target_curve, &target_curve->P, NULL);
            t_ordQ = reduce_to_2_torsion(target_curve, &target_curve->Q, NULL);
            printf("t_ordP=%d, t_ordQ=%d\n", t_ordP, t_ordQ);
            mont_curve_printf(target_curve);
            printf("target P: ");
            mont_pt_printf(&target_curve->P);
            printf("target Q: ");
            mont_pt_printf(&target_curve->Q);
            if (t_ordP < 0 || t_ordQ < 0) {
                rc = 1;
                goto end;
            }

            j_inv(p, target_curve, &raw_key);
            fp2_get_key(&raw_key, item.key);
            fprintf(stderr, "Target key: %s\n", item.key);
            if (hsearch(item, FIND) == NULL) {
                // Add curve to table and queue.
                item.data = malloc(2*sizeof(int));

                if(hsearch(item, ENTER) == NULL) {
                    fprintf(stderr, "Failed to set table value");
                    rc = 1;
                    goto end;
                }
                q[end++] = target_curve;
            } else {
                mont_curve_clear(p, target_curve);
            }
        }
        // Clean up so only wavefront is in memory.
        //free(found_item->data); 
        mont_curve_clear(p, q[front++]);
    }
    
end:
    for (int i = front; i < end; i++) free(q[i]);
    for (int i = 0; i < 3; i++) mont_pt_clear(p, iso2_points[i]);
    fp2_Clear(p, &raw_key);
    
    return rc;
}

int main() {
    int rc = 0;
    unsigned long long num_edges = 0;
    char hkey[H_KEY_SIZE];
    char edges[NUM_NODES][2*H_KEY_SIZE + 5];
    mont_curve_int_t **q = NULL;
    sike_params_t params = { 0 };
    sike_setup_params(&SIKEp33, &params);

    rc = !hcreate(NUM_NODES);
    if (rc) {
        fprintf(stderr, "Failed to initialize table\n");
        goto end;
    }

    q = malloc(NUM_NODES*sizeof(mont_curve_int_t*));
    
    mont_curve_int_t *initial_curve = &params.EA;
    ff_Params *p = initial_curve->ffData;
    q[0] = initial_curve;

    fp2 raw_key = { 0 };
    fp2_Init(p, &raw_key);
    j_inv(p, initial_curve, &raw_key);
    fp2_get_key(&raw_key, hkey);
    int initial_orders[2] = { params.eA, params.eA - 1 };
    ENTRY item = {
        .key = hkey,
        .data = initial_orders,
    };
    hsearch(item, ENTER);

    fprintf(stderr, "Search initialized\n");

    isogeny_bfs(p, q, 1, (char**)edges, &num_edges);
    
end:
    sike_teardown_params(&params);
    hdestroy();
    free(q);
    return rc;
}
