#include <stdio.h>

// external declaration/prototype to call from fortran
void c_wrap_longitude(double *lam, int *radians);

// external declaration/prototype to call from fortran
double c_wrap_latitude(double *phi, int *radians);

int main() {
        int radians = 0;
        double lam = -185.0;
        c_wrap_longitude(&lam, &radians);
        printf("-185 -> %f\n", lam);
        double phi = -135.0;
        double res2 = c_wrap_latitude(&phi, &radians);
        printf("-135 -> %f\n", res2);
        return 0;
}
