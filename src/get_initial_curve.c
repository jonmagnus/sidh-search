#include <sike_params.h>
#include <montgomery.h>
#include <stdlib.h>
#include <stdio.h>
#include <mont_utils.h>

//TODO: Write method for finding linear independant basis with 2^(eA - 1) Q_A = (0,0).
//TODO: Write method to avoid P,Q lying over same 2-torison -- leads to missing edges.

void mount_generic_bases(sike_params_raw_t *partial_raw_params) {
    ff_Params *p = malloc(sizeof(ff_Params));
    set_gmp_fp_params(p);
    fp_ImportHex(partial_raw_params->p, p->mod);
    mont_curve_int_t curve = { 0 };
    mont_curve_init(p, &curve);
    fp2_Set(p, &curve.a, partial_raw_params->A, 0);
    fp2_Set(p, &curve.b, partial_raw_params->B, 0);
    long eA, eB;
    eA = strtol(partial_raw_params->eA, NULL, 0);
    eB = strtol(partial_raw_params->eB, NULL, 0);

    mont_pt_t P = { 0 };
    mont_pt_init(p, &P);
    const int is_alice[4] = { 1, 1, 0, 0 };
    const struct {
        char **x0;
        char **x1;
        char **y0;
        char **y1;
    } pt_hex[] = {
        {
            .x0 = (char**)&partial_raw_params->xQA0,
            .x1 = (char**)&partial_raw_params->xQA1,
            .y0 = (char**)&partial_raw_params->yQA0,
            .y1 = (char**)&partial_raw_params->yQA1,
        },
        {
            .x0 = (char**)&partial_raw_params->xPA0,
            .x1 = (char**)&partial_raw_params->xPA1,
            .y0 = (char**)&partial_raw_params->yPA0,
            .y1 = (char**)&partial_raw_params->yPA1,
        },
        {
            .x0 = (char**)&partial_raw_params->xQB0,
            .x1 = (char**)&partial_raw_params->xQB1,
            .y0 = (char**)&partial_raw_params->yQB0,
            .y1 = (char**)&partial_raw_params->yQB1,
        },
        {
            .x0 = (char**)&partial_raw_params->xPB0,
            .x1 = (char**)&partial_raw_params->xPB1,
            .y0 = (char**)&partial_raw_params->yPB0,
            .y1 = (char**)&partial_raw_params->yPB1,
        },
    };
    for (int i = 0; i < 4; i++) {
        find_basis(&curve, eA, eB, is_alice[i], &P);
        *pt_hex[i].x0 = malloc(19);
        *pt_hex[i].x1 = malloc(19);
        *pt_hex[i].y0 = malloc(19);
        *pt_hex[i].y1 = malloc(19);
        gmp_sprintf(*pt_hex[i].x0, "0x%ZX", P.x.x0);
        gmp_sprintf(*pt_hex[i].x1, "0x%ZX", P.x.x1);
        gmp_sprintf(*pt_hex[i].y0, "0x%ZX", P.y.x0);
        gmp_sprintf(*pt_hex[i].y1, "0x%ZX", P.y.x1);
    }
    mont_pt_clear(p, &P);
    mont_curve_clear(p, &curve);
}

void get_initial_curve(int eA, int eB, sike_params_raw_t *raw_params) {
    unsigned long long p = 1llu << eA;
    for (int i = 0; i < eB; i++) {
        p = (p << 1) + p;
    }
    p--;

    raw_params->p = malloc(19);
    sprintf((char*)raw_params->p, "0x%llX", p);
    raw_params->A = 6;
    raw_params->B = 1;
    raw_params->lA = "2";
    raw_params->lB = "3";
    raw_params->eA = malloc(5);
    raw_params->eB = malloc(5);
    sprintf((char*)raw_params->eA, "%d", eA);
    sprintf((char*)raw_params->eB, "%d", eB);
}
