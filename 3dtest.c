#include "grid/octree.h"
#include "embed.h"
#include "mycompressible_nofr.h"
#include "view.h"
#include "cmap.h"
#include "ibl.h"


#define maxlevel 6
#define minlevel 4

#define mach 0.6
#define rho0 1.
#define p0 1.

int main() {
    L0 = 1.0;
    size(L0);

    init_grid(1 << minlevel);
    origin(-L0/2., -L0/2., -L0/2.);

    run();
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

// 固体边界 - 滑移壁（无穿透）
// rho[embed] = neumann(0.);
// E[embed] = neumann(0.);
// w.n[embed] = dirichlet(0.);  // 无穿透
// w.t[embed] = neumann(0.);    // 滑移

scalar ue[];
scalar w_grad_scalar[];

event solid_init(i = 0) {

    do {
        solid(cs, fs, sq(x) + sq(y) + sq(z) - sq(0.25));
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);
}

event init(i = 0) {
    foreach() {
        if (cs[] > 0.) {
            rho[] = rho0;
            w.x[] = mach * rho0;
            w.y[] = 0.;
            w.z[] = 0.;
            E[] = p0/(gammao - 1.) + 0.5*rho0*sq(mach);
        } else {
            rho[] = 0.;
            w.x[] = 0.;
            w.y[] = 0.;
            w.z[] = 0.;
            E[] = 0.;
        }
    }
}

event adapt(i++) {
    adapt_wavelet({cs, rho, w}, (double[]){1e-30, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2}, maxlevel, minlevel);
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

double scalar_abs_max(scalar s) {
    double max = 0.;
    foreach(reduction(max:max)) {
        if (max < fabs(s[])) {
            max = fabs(s[]);
        }
    }
    return max;
}

event logfile(i++) {
    w_max = scalar_max(w.x);
    cs_min = HUGE;
    foreach(reduction(min:cs_min)) {
        if (cs[] > 0. && cs[] < cs_min) {
            cs_min = cs[];
        }
    }
    fprintf(stderr, "i: %04d t: %.4e dt: %.4e grid: %ld w_max: %.4e cs_min: %.4e\n", i, t, dt, grid->n, w_max, cs_min);
}




event output(t += 1.0) {
    static int frame = 0;
    w_max = scalar_max(w.x);
    w_min = scalar_min(w.x);

    w_min = fmax(0., w_min);

    char filename[80];

    sprintf(filename, "output/grid_%04d.png", frame);
    view(width = 2000, height = 2000);
    clear();
    // squares("w.x", map=viridis, min=w_min, max=w_max, linear=false);
    draw_vof("cs", "fs", color = "w.x");
    colorbar(map=viridis, min=w_min, max=w_max, label="w.x", pos={-0.95, 0.});
    save(filename);
    frame++;
}

event test(t = 1.0) {

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
    cutcell_tangential_gradient_extend(ue, w_grad);

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
    draw_vof("cs", "fs", color = "w_grad.x", map = blue_white_red, min = -scalar_max(w_grad.x), max = scalar_max(w_grad.x), linear = false);
    // squares("w_grad.x", map=blue_white_red, min = -scalar_max(w_grad.x), max = scalar_max(w_grad.x), linear = false);
    colorbar(map = blue_white_red, min = -scalar_max(w_grad.x), max = scalar_max(w_grad.x));
    save("output/w_grad_scalar.png");
}

event end(t = 5.0);