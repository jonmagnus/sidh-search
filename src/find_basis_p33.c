#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <fp.h>
#include <fp2.h>
#include <sike_params.h>
#include <assert.h>

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

void mont_eval_f(const mont_curve_int_t* curve, const fp2* xP, fp2* fxP) {
    const ff_Params *p = curve->ffData;
    const fp2 *a = &curve->a;

    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);
    fp2 *t1 = &T.x, *t2 = &T.y;

    fp2_Square(p, xP, t1);        // t1 = xP^2
    fp2_Multiply(p, xP, t1, t2);  // t2 = xP^3
    fp2_Multiply(p, a, t1, t1);   // t1 = a*xP^2
    fp2_Add(p, t2, t1, t1);       // t1 = xP^3+a*xP^2
    fp2_Add(p, t1, xP, t1);       // t1 = xP^3+a*xP^2+xP
    
    fp2_Copy(p, t1, fxP);

    mont_pt_clear(p, &T);
}

void get_yP(const mont_curve_int_t* curve,
            const fp2* xP,
            mont_pt_t* P) {
    const ff_Params *p = curve->ffData;

    fp2 *t1 = malloc(sizeof(fp2));
    fp2_Init(p, t1);

    mont_eval_f(curve, xP, t1);   // t1 = xP^3+a*xP^2+xP
    fp2_Sqrt(p, t1, t1, 0);       // t1 = sqrt(xP^3+a*xP^2+xP)

    fp2_Copy(p, t1, &P->y);
    fp2_Copy(p, xP, &P->x);

    fp2_Clear(p, t1);
}

int mont_is_on_curve(const mont_curve_int_t* curve, const mont_pt_t *P) {
    const ff_Params *p = curve->ffData;

    const fp2 *xP = &P->x, *yP = &P->y;
    const fp2 *b = &curve->b;

    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);
    fp2 *t1 = &T.x, *t2 = &T.y;

    mont_eval_f(curve, xP, t1);   // t1 = xP^3+a*xP^2+xP
    
    fp2_Square(p, yP, t2);        // t2 = yP^2
    fp2_Multiply(p, t2, b, t2);   // t2 = b*yP^2

    int ans = fp2_IsEqual(p, t1, t2);
    
    mont_pt_clear(p, &T);
    
    return ans;
}

/* returns 1 for True, 0 for False */
static int mont_is_inf_affine(const mont_curve_int_t* curve, const mont_pt_t *P) {
  return (fp2_IsConst( curve->ffData, &P->x, 0, 0)
          && fp2_IsConst( curve->ffData, &P->y, 0, 0 ));
}

int is_square(const ff_Params* p, const fp2* a) {
    fp2* b = malloc(sizeof(fp2));
    fp2_Init(p, b);
    fp2_Sqrt(p, a, b, 0);
    fp2_Square(p, b, b);
    int ans = fp2_IsEqual(p, a, b);
    fp2_Clear(p, b);
    return ans;
}

void mont_rand_pt(const mont_curve_int_t* curve, mont_pt_t* P) {
    const ff_Params *p = curve->ffData;
    const fp2 *a = &curve->a;
    fp2 u = { 0 };
    fp2_Init(p, &u);
    fp2_Set(p, &u, 5134, 1);
    //fp2_Set(p, &u, 4, 1);

    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);
    fp2 *t1 = &T.x, *t2 = &T.y;
    fp2 *t0 = malloc(sizeof(fp2));
    fp2_Init(p, t0);
    fp2_Rand(p, t0);                // t0 = r_ random value

    do {
        fp2_Add(p, t1, t0, t1);         // t1 = r = k*r_ random value

        fp2_Square(p, t1, t1);          // t1 = r^2
        fp2_Multiply(p, &u, t1, t1);    // t1 = u*r^2
        fp2_Set(p, t2, 1, 0);           // t2 = 1
        fp2_Add(p, t1, t2, t1);         // t1 = 1 + u*r^2
        fp2_Invert(p, t1, t1);          // t1 = 1 / (1 + u*r^2)
        fp2_Negative(p, t1, t1);        // t1 = - 1 / (1 + u*r^2)
        fp2_Multiply(p, a, t1, t1);     // t1 = - a / (1 + u*r^2)

        mont_eval_f(curve, t1, t2);
        if (!is_square(p, t2)) {
            fp2_Add(p, a, t1, t1);      // t1 = a - a/(1 + u*r^2)
            fp2_Negative(p, t1, t1);    // t1 = a/(1 + u*r^2) - a
        }
    } while (is_square(p, t1));

    get_yP(curve, t1, P);

    mont_pt_clear(p, &T);
    fp2_Clear(p, &u);
}

void find_basis(const mont_curve_int_t* curve, int eA, int eB, int is_alice, mont_pt_t *P) {
    const ff_Params *p = curve->ffData;
    fp2 xP = { 0 };
    fp2_Init(p, &xP);

    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);

    int order = 0, runs = 0;
    do { 
        mont_rand_pt(curve, P);

        if (!mont_is_on_curve(curve, P)) {
            fprintf(stderr, "P is not on curve\n");
            continue;
        }

        order = 0;
        if (is_alice) {
            xTPLe(curve, P, eB, P);
            mont_pt_copy(p, P, &T);
            while (!mont_is_inf_affine(curve, &T)) {
                xDBL(curve, &T, &T);
                order++;
            }
        } else {
            xDBLe(curve, P, eA, P);
            mont_pt_copy(p, P, &T);
            while (!mont_is_inf_affine(curve, &T)) {
                xTPL(curve, &T, &T);
                order++;
            }
        }
        ++runs;
    } while (order < (is_alice ? eA : eB) && runs < 10);

    mont_pt_clear(p, &T);
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
