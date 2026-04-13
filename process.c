#include <stdio.h>
#include <math.h>

double x[1000], y[1000];

#define AOA (-2.31)

int main() {
    FILE* fp = fopen("airfoil/RAE2822.gnu", "r");
    int n = 0;
    while (fscanf(fp, "%lf %lf", &x[n], &y[n]) != EOF) {
        n++;
    }
    fclose(fp);
    fprintf(stderr, "n: %d\n", n);

    // double alpha = AOA * M_PI / 180.;

    // for (int i = 0; i < n; i++) {
    //     x[i] = x[i] * cos(alpha) - y[i] * sin(alpha) - 0.5;
    //     y[i] = x[i] * sin(alpha) + y[i] * cos(alpha);
    // }

    FILE* fp_out = fopen("airfoil/RAE2822_processed_pointwise.dat", "w");

    for (int i = 64; i >= 0; i--) {
        fprintf(fp_out, "%.6f %.6f %.6f\n", x[i], y[i], 0.);
    }

    for (int i = 66; i < 130; i++) {
        fprintf(fp_out, "%.6f %.6f %.6f\n", x[i], y[i], 0.);
    }
    fclose(fp_out);

    return 0;
}