#include <montgomery.h>
#include <sike_params_small.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

static int mont_is_inf_affine(const mont_curve_int_t* curve, const mont_pt_t *P) {
  return (fp2_IsConst( curve->ffData, &P->x, 0, 0)
          && fp2_IsConst( curve->ffData, &P->y, 0, 0 ));
}

unsigned int find_order(const mont_curve_int_t* curve, int base, const mont_pt_t *P) {
    const ff_Params *p = curve->ffData;
    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);
    mont_pt_copy(p, P, &T);
    unsigned int order = 0;
    if (base == 2) {
        do {
            order++;
            xDBL(curve, &T, &T);
        } while(!mont_is_inf_affine(curve, &T));
    } else if (base == 3) {
        do {
            order++;
            xTPL(curve, &T, &T);
        } while(!mont_is_inf_affine(curve, &T));
    } else order = -1;
    
    mont_pt_clear(p, &T);

    return order;
}

int main() {
    const sike_params_raw_t *raw_params[] = {
        &SIKEp33,
        &SIKEp434,
        &SIKEp503,
        &SIKEp610,
        &SIKEp751};
    for (int i = 0; i < 4; i++) {
        sike_params_t *params = malloc(sizeof(sike_params_t));
        sike_setup_params(raw_params[i],  params);

        const mont_curve_int_t *curveA = &params->EA;
        const mont_curve_int_t *curveB = &params->EB;
        const ff_Params *p = curveA->ffData;

        int orderPA, orderQA;
        orderPA = find_order(curveA, 2, &curveA->P);
        orderQA = find_order(curveA, 2, &curveA->Q);

        if (!i) {
            mont_pt_t dblQ = { 0 };
            mont_pt_init(p, &dblQ);
            xDBL(curveA, &curveA->Q, &dblQ);
            gmp_printf("%Zx\n", dblQ.x.x0);
            gmp_printf("%Zx\n", dblQ.x.x1);
            gmp_printf("%Zx\n", dblQ.y.x0);
            gmp_printf("%Zx\n", dblQ.y.x1);
            mont_pt_clear(p, &dblQ);
        }

        fprintf(stderr, "eA=%lu, orderPA=%d, orderQA=%d\n", params->eA, orderPA, orderQA);

        int orderPB, orderQB;
        orderPB = find_order(curveB, 3, &curveB->P);
        orderQB = find_order(curveB, 3, &curveB->Q);
        
        fprintf(stderr, "eB=%lu, orderPB=%d, orderQB=%d\n", params->eB, orderPB, orderQB);

        sike_teardown_params(params);
    }
}
