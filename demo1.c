#include "grid/quadtree.h"
#include "embed.h"
#include "all-mach.h"
#include "distance.h"
#include "view.h"
#include "cmap.h"

#define maxlevel 10
#define minlevel 6

scalar d[];
coord * points = NULL;

#define mach 0.6
#define rho0 1.
#define p0 1.
#define gammao 1.4

face vector muv[];
scalar rhoc2v[], rhov[];
face vector alphav[];

int main() {
    L0 = 5.0;
    size(L0);
    origin(-L0/2., -L0/2.);

    init_grid(1 << minlevel);
    
    mu = muv;

    run();
}

q.n[left] = dirichlet(mach*rho0);
p[left] = neumann(0.);

q.n[right] = neumann(0.);
p[right] = dirichlet(0.);

q.n[embed] = dirichlet(0.);
q.t[embed] = dirichlet(0.);
p[embed] = neumann(0.);

event defaults(i = 0) {
    rhoc2 = rhoc2v;
    rho = rhov;
    alpha = alphav;
}

event init_solid(i = 0) {
    FILE* fp = fopen("NACA0012_blunt.gnu", "r");
    points = input_xy(fp);
    distance(d, points);
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);

    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);
}

event properties(i++) {
    foreach_face() {
        muv.x[] = fm.x[] * 1e-3;
        alphav.x[] = fm.x[]/rho0;
    }
    foreach() {
        rhov[] = rho0 * cm[];
        ps[] = p[];
        rhoc2v[] = gammao * p[];
    }
}

event adapt(i++) {
    adapt_wavelet({cs, q}, (double[]){1e-30, 1e-2, 1e-2}, maxlevel, minlevel);
}

event logfile(i++) {
    fprintf(stderr, "i: %04d t: %.4e dt: %.4e grid: %ld\n", i, t, dt, grid->n);
}

event output(t += 1.0) {
    char name[80];
    static int count = 0;
    sprintf(name, "output/grid_%04d.png", count++);
    view(width = 2000, height = 2000);
    clear();
    squares("q.x", min=0., max=1.5, map=viridis, linear=false);
    save(name);
}

event end(t = 5.0);