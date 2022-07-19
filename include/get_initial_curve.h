#ifndef GET_INITIAL_CURVE_H
#define GET_INITIAL_CURVE_H

#include <sike_params.h>

void mount_generic_bases(sike_params_raw_t *raw_params);

void get_initial_curve(int eA, int eB, sike_params_raw_t *raw_params);

#endif // GET_INITIAL_CURVE_H
