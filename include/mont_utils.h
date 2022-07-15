#ifndef MONT_UTILS_H
#define MONT_UTILS_H

#include <montgomery.h>

void mont_eval_f(const mont_curve_int_t *curve, const fp2 *xP, fp2 *fxP);

int mont_is_on_curve(const mont_curve_int_t *curve, mont_pt_t *P);

int is_square(const ff_Params *p, const fp2 *a);

void mont_get_yP(const mont_curve_int_t *curve,
                 const fp2 *xP,
                 mont_pt_t *P);

void mont_rand_pt(const mont_curve_int_t *curve, mont_pt_t *P);

int mont_is_inf(const ff_Params *p, const mont_pt_t *P);

int mont_is_principal_2_torsion(const mont_curve_int_t *currve,
                                const mont_pt_t *P);

int reduce_to_2_torsion(const mont_curve_int_t *curve,
                        const mont_pt_t *P,
                        mont_pt_t *V);

#endif // MONT_UTILS_H
