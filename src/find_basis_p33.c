#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <fp.h>
#include <fp2.h>
#include <sike_params.h>
#include <assert.h>
#include <mont_utils.h>

const struct {
    const char* p;
    
    int eA;
    int eB;
} p33 = {
    .p = "0x01037e1fff",
    .eA = 13,
    .eB = 12,
};

char* fp2_get_str(char* str, int base, const fp2* a) {
    char* x0_str = mpz_get_str(NULL, base, a->x0);
    char* x1_str = mpz_get_str(NULL, base, a->x1);
    if (str == NULL) {
        str = malloc((strlen(x0_str) + strlen(x1_str) + 5)*sizeof(char));
    }
    sprintf(str, "(%s, %s)", x0_str, x1_str);
    return str;
}

char* mont_pt_get_str(char* str, int base, const mont_pt_t* P) {
    char* x_str = fp2_get_str(NULL, base, &P->x);
    char* y_str = fp2_get_str(NULL, base, &P->y);
    if (str == NULL) {
        str = malloc((strlen(x_str) + strlen(y_str) + 2)*sizeof(char));
    }
    sprintf(str, "%s %s", x_str, y_str);
    return str;
}

int main() {

    ff_Params* p = malloc(sizeof(ff_Params));
    set_gmp_fp_params(p);
   
    fp_Init(p, p->mod);
    fp_ImportHex(p33.p, p->mod);

    fp2 A = { 0 }, B = { 0 };
    fp2_Init(p, &A);
    fp2_Init(p, &B);
    fp2_Set(p, &A, 6, 0);       // Use same curve as every other spec
    fp2_Set(p, &B, 1, 0);       // Isomorphic montgomery curve

    /////////////////////
    // Find basis points.
    /////////////////////
    printf(
            "A=%s\nB=%s, square=%d\n",
            mpz_get_str(NULL, 0, A.x0),
            mpz_get_str(NULL, 0, B.x0),
            is_square(p, &B));

    mont_curve_int_t curve = {
        .ffData = p,
        .a = A,
        .b = B,
    };

    mont_pt_t P = { 0 };
    mont_pt_init(p, &P);
    find_basis(&curve, p33.eA, p33.eB, 1, &P);
    printf("PA=%s\n", mont_pt_get_str(NULL, 16, &P));
    find_basis(&curve, p33.eA, p33.eB, 1, &P);
    xDBL(&curve, &P, &P);       // Q needs to be a double.
    printf("QA=%s\n", mont_pt_get_str(NULL, 16, &P));
    find_basis(&curve, p33.eA, p33.eB, 0, &P);
    printf("PB=%s\n", mont_pt_get_str(NULL, 16, &P));
    find_basis(&curve, p33.eA, p33.eB, 0, &P);
    printf("QB=%s\n", mont_pt_get_str(NULL, 16, &P));

    fp_Clear(p, p->mod);
    fp2_Clear(p, &A);
    fp2_Clear(p, &B);
    
    free(p);
    return 0;
};
