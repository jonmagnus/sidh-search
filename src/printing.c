#include <gmp.h>
#include <montgomery.h>

void fp2_printf(const fp2 *a) {
    gmp_printf("%Zu + i*%Zu\n",
               a->x0,
               a->x1);
}

void mont_pt_printf(const mont_pt_t *P) {
    gmp_printf("(%Zu + i*%Zu, %Zu + i*%Zu)\n",
               P->x.x0,
               P->x.x1,
               P->y.x0,
               P->y.x1);
}

void mont_curve_printf(const mont_curve_int_t *curve) {
    gmp_printf("a: %Zu + i*%Zu\n",
               curve->a.x0,
               curve->a.x1);
    gmp_printf("b: %Zu + i*%Zu\n",
               curve->b.x0,
               curve->b.x1);
}

