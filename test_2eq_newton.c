#include "grid/quadtree.h"
#include "embed.h"
#include "distance.h"
#include "mycompressible_nofr.h"
#include "view.h"
#include "cmap.h"
#include "ibl.h"

#define maxlevel 10
#define minlevel 6

#define mach 0.6
#define rho0 1.
#define p0 (1./1.4)

scalar d[];
coord* p = NULL;

#define ibl_niter 1000000
#define output_interval 1000

scalar ibl_theta[];
scalar ibl_H[];
scalar Cf2[];
scalar delta[];

scalar ue[];
vector ue_vec[];
scalar Me[];

double stag_x, stag_y, ue_min = HUGE, tail_x = -HUGE;

void update_ue();
void output_data();
void ibl_solver();

int main() {
    L0 = 5.0;
    size(L0);

    init_grid(1 << minlevel);
    origin(-L0/2., -L0/2.);

    CFL = 0.3;

    restore("restart");
    boundary(all);

    FILE* fp = fopen("NACA0012_blunt_new.gnu", "r");
    p = input_xy(fp);
    distance(d, p);
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);

    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    build_cutcell_cache();
    
    foreach_cache(cutcells) {
        if (ue[] < ue_min) {
            ue_min = ue[];
            stag_x = x;
            stag_y = y;
        }
        ibl_theta[] = 1e-6;
        ibl_H[] = 1.4;
        if (x > tail_x) {
            tail_x = x;
        }
    }
    update_ue();

    for (int iter = 0; iter < ibl_niter; iter++) {

        ibl_solver();

        if (iter % 100 == 0) 
            fprintf(stderr, "[log] iter:%04d \n", iter);

        if (iter % output_interval == 0)        
            output_data();
            
    }
}

#include "ibl_solver.h"

void update_ue() {
    foreach_cache(cutcells) {
        coord n, b;
        embed_geometry(point, &b, &n);

        double rho_s = my_embed_interpolate(point, rho, b);
        double ux = my_embed_interpolate(point, w.x, b) / rho_s;
        double uy = my_embed_interpolate(point, w.y, b) / rho_s;
        double un = ux * n.x + uy * n.y;
        double etx = ux - un * n.x;
        double ety = uy - un * n.y;
        ue[] = sqrt(etx * etx + ety * ety);
        double norm = sqrt(sq(etx) + sq(ety));
        ue_vec.x[] = etx / norm;
        ue_vec.y[] = ety / norm;

        double E_s = my_embed_interpolate(point, E, b);
        double p_s = (gammao - 1.) * (E_s - 0.5*(sq(ux) + sq(uy))/rho_s);
        double a_s = sqrt(gammao * p_s / rho_s);
        Me[] = ue[] / fmax(a_s, 1e-10);
    }
    boundary({ue, ue_vec, Me});
}

void output_data() {
    FILE* f_theta = fopen("data/theta.dat", "w");
    fprintf(f_theta, "# x theta\n");
    foreach_cache(cutcells) {
        if (y > 0.)
            fprintf(f_theta, "%.6e %.6e\n", (x - stag_x)/0.97, ibl_theta[]);
    }
    fflush(f_theta);
    fclose(f_theta);
    system("gnuplot plot_th.gp");

    FILE* f_Cf = fopen("data/Cf2.dat", "w");
    fprintf(f_Cf, "# x Cf2\n");
    foreach_cache(cutcells) {
        if (y > 0.)
            fprintf(f_Cf, "%.6e %.6e\n", x - stag_x, Cf2[]*2.);
    }
    fflush(f_Cf);
    fclose(f_Cf);
    system("gnuplot plot_Cf.gp");

    FILE* f_H = fopen("data/H.dat", "w");
    fprintf(f_H, "# x H\n");
    foreach_cache(cutcells) {
        if (y > 0.)
            fprintf(f_H, "%.6e %.6e\n", x - stag_x, ibl_H[]);
    }
    fflush(f_H);
    fclose(f_H);
    system("gnuplot plot_H.gp");

    FILE* f_ue = fopen("data/ue.dat", "w");
    fprintf(f_ue, "# x ue\n");
    foreach_cache(cutcells) {
        if (y > 0.)
            fprintf(f_ue, "%.6e %.6e\n", x - stag_x, ue[]/0.6);
    }
    fflush(f_ue);
    fclose(f_ue);
    system("gnuplot plot_ue.gp");

    FILE* f_delta = fopen("data/delta.dat", "w");
    fprintf(f_delta, "# x delta\n");
    foreach_cache(cutcells) {
        if (y > 0.)
            fprintf(f_delta, "%.6e %.6e\n", x - stag_x, delta[]);
    }
    fflush(f_delta);
    fclose(f_delta);
    system("gnuplot plot_delta.gp");
}