#include <stdio.h>
#include <stdlib.h>

#define x_offset -1.0

int main() {
    FILE * fp = fopen("NACA0012.gnu", "r");
    int count = 0;
    double x, y;
    while (fscanf(fp, "%lf %lf", &x, &y) != EOF) {
        count++;
    }

    printf("count: %d\n", count);
    double* x_list = (double*) malloc(count * sizeof(double));
    double* y_list = (double*) malloc(count * sizeof(double));

    rewind(fp);

    fscanf(fp, "%lf %lf", &x_list[0], &y_list[0]);
    for (int i = 1; i < count; i++) {
        fscanf(fp, "%lf %lf", &x_list[i], &y_list[i]);
    }

    FILE * fp2 = fopen("NACA0012_processed.gnu", "w");
    for (int i = 0; i < count; i++) {
        fprintf(fp2, "%lf %lf\n", x_list[i] + x_offset, y_list[i]);
    }
    fflush(fp2);
    fclose(fp2);
    free(x_list);
    free(y_list);
    fclose(fp);
}