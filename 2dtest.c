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
#define p0 1.

scalar d[];
coord* p = NULL;

scalar * my_dumplist = {rho, w.x, w.y, E, cs};

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
// rho[embed] = neumann(0.);
// E[embed] = neumann(0.);
// w.n[embed] = dirichlet(0.);  // 无穿透
// w.t[embed] = neumann(0.);    // 滑移

scalar ue[];

event solid_init(i = 0) {
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

event logfile(i++) {
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
    fprintf(stderr, "i: %04d t: %.4e dt: %.4e grid: %ld w_max: %.4e cs_min: %.4e dw: %.4e\n", i, t, dt, grid->n, w_max, cs_min, dw);
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

event end(t = 5.0) {
    if (p) {
        free (p);
    }

    d.delete = NULL;

    dump("restart", my_dumplist);
}
