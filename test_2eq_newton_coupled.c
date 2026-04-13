#include "grid/quadtree.h"
#include "embed.h"
#include "distance.h"
#include "mycompressible_nofr.h"
#include "view.h"
#include "cmap.h"
#include "ibl.h"

#define maxlevel 10
#define minlevel 6

#define mach 0.729
#define rho0 1.
#define p0 1.
#define ibl_Re (6.5e6)

#define ibl_nu (rho0 * u_inf / ibl_Re)

#define u_inf (mach * sqrt(gammao * p0 / rho0))
#define E_inf (p0/(gammao - 1.) + 0.5 * rho0 * sq(u_inf))
#define w_inf (rho0 * u_inf)

scalar d[];
coord* p = NULL;

#define ibl_niter 10
#define output_interval 100

scalar ibl_theta[];
scalar ibl_H[];
scalar Cf2[];
scalar delta_star[];
scalar ibl_Ctau[];

scalar rho_u_delta[];

scalar ue[];
vector ue_vec[];
scalar Me[];

double stag_x = HUGE, stag_y, ue_min = HUGE, tail_x = -HUGE;

void update_ue();
void output_data();
void ibl_solver();

int main() {
    L0 = 5.0;
    size(L0);

    init_grid(1 << minlevel);
    origin(-L0/2., -L0/2.);

    CFL = 0.3;

    run();
}

event init_solid(i = 0) {
    
    // FILE* fp = fopen("airfoil/NACA0012_blunt_new.gnu", "r");
    FILE* fp = fopen("airfoil/RAE2822_processed.gnu", "r");
    p = input_xy(fp);
    distance(d, p);
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);

    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    build_cutcell_cache();
    
    foreach() {
        if (cs[] > 0.) {
            rho[] = rho0;
            w.x[] = w_inf;
            w.y[] = 0.;
            E[] = E_inf;
        } else {
            rho[] = 0.;
            w.x[] = 0.;
            w.y[] = 0.;
            E[] = 0.;
        }
    }

    foreach_cache(cutcells) {
        if (x < stag_x) {
            stag_x = x;
            stag_y = y;
        }
        ibl_theta[] = 1e-6;
        ibl_H[] = 2.2; 
        ibl_Ctau[] = 0.005;
        if (x > tail_x) {
            tail_x = x;
        }
    }
}

rho[left] = dirichlet(rho0);
E[left] = dirichlet(E_inf);
w.x[left] = dirichlet(w_inf);
w.y[left] = dirichlet(0.);

// 流出边界
rho[right] = neumann(0.);
E[right] = neumann(0.);
w.x[right] = neumann(0.);
w.y[right] = neumann(0.);

event adapt(i++) {
    adapt_wavelet({cs, rho, w}, (double[]){1e-30, 1e-2, 1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
    fractions_cleanup(cs, fs);
    build_cutcell_cache();
}


event end(t = 10.);

#include "ibl_solver.h"

event logfile(i += 10) {
    fprintf(stderr, "[log] i: %04d t: %.4e dt: %.4e grid: %ld\n", i, t, dt, grid->n);
}

event ibl_solver(i++) {
    update_ue();
    for (int iter = 0; iter < ibl_niter; iter++) {
        ibl_solver();
    }

    if (i % output_interval == 0)
        output_data();

    // foreach_cache(cutcells) 
    //     rho_u_delta[] = rho[] * ue[] * delta_star[];
    
    // if (t > 0.1) {
    //     foreach_cache(cutcells) {
    //         if (fabs(x - tail_x) < 0.02) continue;
    //         coord n, b;
    //         coord e = {ue_vec.x[], ue_vec.y[]};
    //         double area = embed_geometry(point, &b, &n);
    //         double m = upwind_grad(rho_u_delta, e, point);

    //         double w2 = 0.;
    //         foreach_dimension()
    //             w2 += sq(w.x[]);

    //         double p_e = (gammao - 1.) * (E[] - 0.5 * w2 / rho[]);
    //         double H_e = (E[] + p_e) / rho[];

    //         double factor = m * area / (fmax(cs[], CS_FLOOR) * Delta);

    //         rho[] += factor * dt;
    //         foreach_dimension() 
    //             w.x[] += factor * ue_vec.x[] * ue[] * dt;
            
    //         E[] += factor * H_e * dt;
    //     }
    // }
}

event movie(t += 0.5) {
    double rho_min = 1e30, rho_max = -1e30;
    foreach(reduction(min:rho_min) reduction(max:rho_max)) {
        if (!is_solid) {
            if (rho[] < rho_min) rho_min = rho[];
            if (rho[] > rho_max) rho_max = rho[];
        }
    }
    static int count = 0;
    char name[80];
    sprintf(name, "output/movie_%04d.png", count++);
    view(width = 2000, height = 2000);
    clear();
    squares("rho", map = viridis, min = rho_min, max = rho_max, linear = false);
    draw_vof("cs", "fs", lc = {0., 0., 0.}, lw = 1.0);
    colorbar(map = viridis, min = rho_min, max = rho_max, label = "rho", pos = {-0.95, 0.});
    save(name);
}

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
        // if (y > 0.) 
            fprintf(f_theta, "%.6e %.6e\n", (x - stag_x)/0.97, ibl_theta[]);
    }
    fflush(f_theta);
    fclose(f_theta);
    system("gnuplot plot_th.gp");

    FILE* f_Cf = fopen("data/Cf2.dat", "w");
    fprintf(f_Cf, "# x Cf2\n");
    foreach_cache(cutcells) {
        // if (y > 0.) 
            fprintf(f_Cf, "%.6e %.6e\n", x - stag_x, 2. * Cf2[]);
    }
    fflush(f_Cf);
    fclose(f_Cf);
    system("gnuplot plot_Cf.gp");

    FILE* f_H = fopen("data/H.dat", "w");
    fprintf(f_H, "# x H\n");
    foreach_cache(cutcells) {
        // if (y > 0.)
            fprintf(f_H, "%.6e %.6e\n", x - stag_x, Hk[]);
    }
    fflush(f_H);
    fclose(f_H);
    system("gnuplot plot_H.gp");

    FILE* f_ue = fopen("data/ue.dat", "w");
    fprintf(f_ue, "# x ue\n");
    foreach_cache(cutcells) {
        // if (y > 0.)
            fprintf(f_ue, "%.6e %.6e\n", x - stag_x, ue[]/u_inf);
    }
    fflush(f_ue);
    fclose(f_ue);
    system("gnuplot plot_ue.gp");

    FILE* f_delta = fopen("data/delta.dat", "w");
    fprintf(f_delta, "# x delta\n");
    foreach_cache(cutcells) {
        // if (y > 0.)
            fprintf(f_delta, "%.6e %.6e\n", x - stag_x, delta_star[]);
    }
    fflush(f_delta);
    fclose(f_delta);
    system("gnuplot plot_delta.gp");

    FILE* f_Ctau = fopen("data/Ctau.dat", "w");
    fprintf(f_Ctau, "# x Ctau\n");
    foreach_cache(cutcells) {
        // if (y > 0.)
            fprintf(f_Ctau, "%.6e %.6e\n", x - stag_x, ibl_Ctau[]);
    }
    fflush(f_Ctau);
    fclose(f_Ctau);
    system("gnuplot plot_Ctau.gp");

    FILE* f_Re_th = fopen("data/Re_th.dat", "w");
    fprintf(f_Re_th, "# x Re_th\n");
    foreach_cache(cutcells) {
        // if (y > 0.)
            fprintf(f_Re_th, "%.6e %.6e\n", x - stag_x, Re_th[]);
    }
    fflush(f_Re_th);
    fclose(f_Re_th);
    system("gnuplot plot_Re_th.gp");
}