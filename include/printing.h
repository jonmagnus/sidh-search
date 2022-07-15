#ifndef PRINTING_H
#define PRINTING_H

#include <montgomery.h>

void fp2_printf(const fp2 *a);

void mont_pt_printf(const mont_pt_t *P);

void mont_curve_printf(const mont_curve_int_t *curve);

#endif // PRINTING_H
