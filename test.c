#include "grid/octree.h"
#include "embed.h"
#include "mycompressible_nofr.h"
#include "view.h"
#include "cmap.h"
#include "ibl.h"

int maxlevel = 6;
int minlevel = 4;

#define mach 0.6
#define rho0 1.
#define p0 1.

double w_max = 0.;
double w_min = 0.;
double rho_max = 0.;
double rho_min = 0.;
double cs_min = 0.;

scalar ue[];
scalar w_grad_scalar[];

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

double scalar_abs_max(scalar s) {
    double max = 0.;
    foreach(reduction(max:max)) {
        if (max < fabs(s[])) {
            max = fabs(s[]);
        }
    }
    return max;
}

int main(int argc, char** argv) {
    restore("restart");
    do {
        solid(cs, fs, sq(x) + sq(y) + sq(z) - sq(0.25));
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    vector w_grad[];
    foreach() {
        foreach_dimension()
            w_grad.x[] = 0.;
        ue[] = 0.;
        w_grad_scalar[] = 0.;
    }

    build_cutcell_cache();
    foreach_cache(cutcells) {
        coord n, b;

        embed_geometry(point, &b, &n);
        double wx = my_embed_interpolate(point, w.x, b);
        double wy = my_embed_interpolate(point, w.y, b);
        double wz = my_embed_interpolate(point, w.z, b);
        double rho_s = my_embed_interpolate(point, rho, b);
        double ux = wx/rho_s, uy = wy/rho_s, uz = wz/rho_s;
        double un = ux*n.x + uy*n.y + uz*n.z;

        ue[] = sqrt(sq(ux - un * n.x) + sq(uy - un * n.y) + sq(uz - un * n.z));
    }
    // cutcell_tangential_gradient(ue, w_grad);
    // cutcell_tangential_gradient_extend(ue, w_grad);
    // cutcell_tangential_gradient_implicit_extend(ue, w_grad);

    // cutcell_tangential_gradient_lsq_forward_fixed(ue, w_grad);
    cutcell_tangential_gradient_lsq_backward(ue, w_grad);

    foreach_cache(cutcells) {
        coord n, b;
        embed_geometry(point, &b, &n);
        foreach_dimension() {
            w_grad_scalar[] += sq(w_grad.x[]);
        }
        w_grad_scalar[] = sqrt(w_grad_scalar[]);
    }

    view(width = 2000, height = 2000);
    clear();
    draw_vof("cs", "fs", color = "w_grad.x", map = blue_white_red, min = -scalar_abs_max(w_grad.x), max = scalar_abs_max(w_grad.x), linear = false);
    colorbar(map = blue_white_red, min = -scalar_abs_max(w_grad.x), max = scalar_abs_max(w_grad.x));
    save("output/w_grad.png");
}

rho[left] = dirichlet(rho0);
E[left] = dirichlet(p0/(gammao - 1.) + 0.5*rho0*sq(mach));
w.x[left] = dirichlet(mach*rho0);
w.y[left] = dirichlet(0.);
w.z[left] = dirichlet(0.);

// 流出边界
rho[right] = neumann(0.);
E[right] = neumann(0.);
w.x[right] = neumann(0.);
w.y[right] = neumann(0.);
w.z[right] = neumann(0.);
