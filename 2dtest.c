#include "grid/quadtree.h"
#include "embed.h"
#include "distance.h"
#include "mycompressible_nofr.h"
// #include "mycompressible_mirror.h"
#include "view.h"
#include "cmap.h"
#include "ibl.h"


#define maxlevel 10
#define minlevel 6

#define mach 0.6
#define rho0 1.
#define p0 1.

scalar d[];
coord* p = NULL;

scalar * my_dumplist = {rho, w.x, w.y, E, cs};

scalar ue[];
vector ue_grad[];

int main() {
    L0 = 5.0;
    size(L0);

    init_grid(1 << minlevel);
    origin(-L0/2., -L0/2.);

    run();
}

rho[left] = dirichlet(rho0);
E[left] = dirichlet(p0/(gammao - 1.) + 0.5*rho0*sq(mach));
w.x[left] = dirichlet(mach*rho0);
w.y[left] = dirichlet(0.);

// 流出边界
rho[right] = neumann(0.);
E[right] = neumann(0.);
w.x[right] = neumann(0.);
w.y[right] = neumann(0.);

// 固体边界 - 滑移壁（无穿透）
rho[embed] = neumann(0.);
E[embed] = neumann(0.);
w.n[embed] = dirichlet(0.);  // 无穿透
w.t[embed] = neumann(0.);    // 滑移

event solid_init(i = 0) {
    // FILE* fp = fopen("NACA0012_blunt.gnu", "r");
    FILE* fp = fopen("NACA0012_processed.gnu", "r");
    p = input_xy(fp);
    distance(d, p);
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);

    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);
}

event init(i = 0) {
    foreach() {
        if (cs[] > 0.) {
            rho[] = rho0;
            w.x[] = mach * rho0;
            w.y[] = 0.;
            E[] = p0/(gammao - 1.) + 0.5*rho0*sq(mach);
        } else {
            rho[] = 0.;
            w.x[] = 0.;
            w.y[] = 0.;
            E[] = 0.;
        }
    }
}

event adapt(i++) {
    adapt_wavelet({cs, rho, w}, (double[]){1e-30, 1e-2, 1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
    fractions_cleanup(cs, fs);
    build_cutcell_cache();
}

double w_max = 0.;
double w_min = 0.;
double rho_max = 0.;
double rho_min = 0.;
double cs_min = 0.;

double scalar_max(scalar s) {
    double max = -HUGE;
    foreach(reduction(max:max)) {
        if (max < s[]) {
            max = s[];
        }
    }
    return max;
}

double scalar_min(scalar s) {
    double min = HUGE;
    foreach(reduction(min:min)) {
        if (min > s[]) {
            min = s[];
        }
    }
    return min;
}

scalar wn[];

event logfile(i += 100) {
    double dw = 0.;
    if (dt > 1e-10)
        dw = change (w.x, wn)/dt;
    else
        dw = 0.;
    w_max = scalar_max(w.x);
    cs_min = HUGE;
    foreach(reduction(min:cs_min)) {
        if (cs[] > 0. && cs[] < cs_min) {
            cs_min = cs[];
        }
    }
    fprintf(stderr, "[log] i: %04d t: %.4e dt: %.4e grid: %ld w_max: %.4e cs_min: %.4e dw: %.4e\n", i, t, dt, grid->n, w_max, cs_min, dw);
}

// event output(t += 1.0) {
//     static int frame = 0;
//     w_max = scalar_max(w.x);
//     w_min = scalar_min(w.x);

//     w_min = fmax(0., w_min);

//     char filename[80];

//     sprintf(filename, "output/grid_%04d.png", frame);
//     view(width = 2000, height = 2000, samples=8);
//     clear();
//     squares("w.x", map=viridis, min=w_min, max=w_max, linear=false);
//     draw_vof("cs", "fs", lc={0., 0., 0.}, lw=1.0);
//     colorbar(map=viridis, min=w_min, max=w_max, label="w.x", pos={-0.95, 0.});
//     save(filename);
//     frame++;
// }

// event test(t = 1.0) {

//     vector w_grad[];
//     foreach() {
//         foreach_dimension()
//             w_grad.x[] = 0.;
//         ue[] = 0.;
//     }

//     build_cutcell_cache();
//     foreach_cache(cutcells) {
//         coord n, b;

//         embed_geometry(point, &b, &n);
//         double wx = embed_interpolate(point, w.x, b);
//         double wy = embed_interpolate(point, w.y, b);
//         double rho_s = embed_interpolate(point, rho, b);
//         double ux = wx/rho_s, uy = wy/rho_s;
//         double un = ux*n.x + uy*n.y;

//         ue[] = sqrt(sq(ux - un * n.x) + sq(uy - un * n.y));
//     }
//     cutcell_tangential_gradient_lsq_backward(ue, w_grad);

//     view();
//     clear();
//     squares("w_grad.x", map=blue_white_red, min = -scalar_max(w_grad.x), max = scalar_max(w_grad.x), linear = false);
//     colorbar(map = blue_white_red, min = -scalar_max(w_grad.x), max = scalar_max(w_grad.x));
//     save("output/w_grad.png");
// }

// event output(t += 0.1) {
//     double x_min = HUGE;
//     foreach_cache(cutcells){
//         if (x < x_min) {
//             x_min = x;
//         }
//     }

//     foreach() {
//         ue[] = 0.;
//         foreach_dimension()
//             ue_grad.x[] = 0.;
//     }

//     foreach_cache(cutcells) {
//         coord n, b;
//         embed_geometry(point, &b, &n);
//         double wx = embed_interpolate(point, w.x, b);
//         double wy = embed_interpolate(point, w.y, b);
//         double rho_s = embed_interpolate(point, rho, b);

//         ue[] = sqrt(sq(wx/rho_s/mach) +sq(wy/rho_s/mach));
//     }

//     FILE * f_ue = fopen("data/ue.dat", "w");
//     fprintf(f_ue, "# x ue\n");
//     foreach_cache(cutcells) {
//         if (y > 0.) {
//             fprintf(f_ue, "%.6e %.6e\n", x - x_min, ue[]);
//         }
//     }
//     fflush(f_ue);
//     fclose(f_ue);

//     system("gnuplot plot_ue.gp");

//     cutcell_tangential_gradient_lsq_backward(ue, ue_grad);
//     // cutcell_tangential_gradient_extend(ue, ue_grad);
//     // cutcell_tangential_gradient(ue, ue_grad);

//     FILE* f_ug = fopen("data/ue_grad.dat", "w");
//     fprintf(f_ug, "# x ue_grad\n");
//     foreach_cache(cutcells) {
//         coord n, b;
//         embed_geometry(point, &b, &n);
//         if (y > 0.) {
//             double value = 0.;
//             coord n, b;
//             embed_geometry(point, &b, &n);
//             foreach_dimension()
//                 value += sq(embed_interpolate(point, ue_grad.x, b));
//             value = sqrt(value);
//             fprintf(f_ug, "%.6e %.6e\n", x - x_min, value);
//         }
//     }
//     fflush(f_ug);
//     fclose(f_ug);

//     system("gnuplot plot_ug.gp");


//     view();
//     clear();
//     draw_vof("cs", "fs", lc={0., 0., 0.}, lw=1.0);
//     save("output/cs.png");
// }

// event blow_wind (i++) {
//     double x_min = HUGE;
//     foreach_cache(cutcells){
//         if (x < x_min) {
//             x_min = x;
//         }
//     }
//     double blow_wind_speed = 0.1;
//     foreach_cache(cutcells) {
//         coord n, b;
//         double area = embed_geometry(point, &b, &n);
//         foreach_dimension()
//             n.x *= -1.0;

//         double coef = 1. + sq(1. + x - x_min);
//         // double coef = 1.0;
//         foreach_dimension()
//             w.x[] += n.x * blow_wind_speed * coef * dt * area / (cm[] * Delta);
//     }
// }

event end(t = 5.0) {
    if (p) {
        free (p);
    }

    d.delete = NULL;

    dump("restart", my_dumplist);
}
