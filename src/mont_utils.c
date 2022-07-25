#include <montgomery.h>
#include <stdlib.h>
#include <stdio.h>
#include <printing.h>

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

int is_square(const ff_Params* p, const fp2* a) {
    fp2* b = malloc(sizeof(fp2));
    fp2_Init(p, b);
    fp2_Sqrt(p, a, b, 0);
    fp2_Square(p, b, b);
    int ans = fp2_IsEqual(p, a, b);
    fp2_Clear(p, b);
    return ans;
}

void mont_get_yP(const mont_curve_int_t* curve,
            const fp2* xP,
            mont_pt_t* P) {
    const ff_Params *p = curve->ffData;

    fp2 *t1, *invb;
    t1 = malloc(sizeof(fp2));
    fp2_Init(p, t1);
    invb = malloc(sizeof(fp2));
    fp2_Init(p, invb);

    fp2_Invert(p, &curve->b, invb);

    mont_eval_f(curve, xP, t1);     // t1 = xP^3+a*xP^2+xP
    fp2_Multiply(p, t1, invb, t1);  // t1 = (xP^3+a*xP^2+xP) / b
    fp2_Sqrt(p, t1, t1, 0);         // t1 = sqrt((xP^3+a*xP^2+xP) / b)

    fp2_Copy(p, t1, &P->y);
    fp2_Copy(p, xP, &P->x);

    fp2_Clear(p, t1);
    fp2_Clear(p, invb);
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

    mont_get_yP(curve, t1, P);

    mont_pt_clear(p, &T);
    fp2_Clear(p, &u);
}

// Might be confused with the point (0:0:1) which is a 2-torsion point.
int mont_is_inf(const ff_Params *p, const mont_pt_t *P) {
    return fp2_IsConst(p, &P->x, 0, 0) && fp2_IsConst(p, &P->y, 0, 0);
}

int mont_is_principal_2_torsion(const mont_curve_int_t *curve,
                                const mont_pt_t *P) {
    const ff_Params *p = curve->ffData;
    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);
    mont_pt_copy(p, P, &T);
    int ans = mont_is_inf(p, &T);
    xDBL(curve, &T, &T);
    ans = !ans && mont_is_inf(p, &T);
    mont_pt_clear(p, &T);
    return ans;
}

int reduce_to_2_torsion(const mont_curve_int_t *curve,
                        const mont_pt_t *P,
                        mont_pt_t *V) {
    const ff_Params *p = curve->ffData;
    mont_pt_t *T;
    T = malloc(sizeof(mont_pt_t));
    mont_pt_init(p, T);
    mont_pt_copy(p, P, T);
    int ord = 1;
    while (!fp2_IsConst(p, &T->y, 0, 0)) {
        xDBL(curve, T, T);
        ord++;
    }

    if (V != NULL) {
        mont_pt_copy(p, T, V);
    }
    mont_pt_clear(p, T);

    return ord;
}

int reduce_to_3_torsion(const mont_curve_int_t *curve,
                        const mont_pt_t *P,
                        mont_pt_t *V) {
    const ff_Params *p = curve->ffData;
    mont_pt_t T = { 0 }, U = { 0 };
    mont_pt_init(p, &T);
    mont_pt_init(p, &U);
    mont_pt_copy(p, P, &T);
    int ord = 0;
    while (!fp2_IsConst(p, &T.y, 0, 0)) {
        mont_pt_copy(p, &T, &U);
        xTPL(curve, &T, &T);
        ord++;
    }

    if (V != NULL) {
        mont_pt_copy(p, &U, V);
    }
    mont_pt_clear(p, &T);
    mont_pt_clear(p, &U);

    return ord;
}

void find_basis(const mont_curve_int_t* curve,
                int eA,
                int eB,
                int is_alice,
                mont_pt_t *P) {
    const ff_Params *p = curve->ffData;

    mont_pt_t T = { 0 };
    mont_pt_init(p, &T);

    int order = 0, runs = 0;
    do { 
        if (runs > 100) exit(1);
        ++runs;
        mont_rand_pt(curve, P);

        if (!mont_is_on_curve(curve, P)) {
            continue;
        }

        order = 0;
        if (is_alice) {
            xTPLe(curve, P, eB, P);
            mont_pt_copy(p, P, &T);
            while (!fp2_IsConst(p, &T.y, 0, 0)) {
                xDBL(curve, &T, &T);
                order++;
            }
            order++;
        } else {
            xDBLe(curve, P, eA, P);
            mont_pt_copy(p, P, &T);
            while (!fp2_IsConst(p, &T.x, 0, 0)) {
                xTPL(curve, &T, &T);
                order++;
            }
        }
    } while (order < (is_alice ? eA : eB));

    mont_pt_clear(p, &T);
}

