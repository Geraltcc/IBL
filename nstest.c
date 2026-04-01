#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "distance.h"
#include "view.h"
#include "cmap.h"

#define maxlevel 10
#define minlevel 6
#define U 1.0
#define MU 1e-3

scalar d[];
coord* points = NULL;

face vector muv[];

int main() {
    L0 = 5.0;
    size(L0);
    origin(-L0/2., -L0/2.);

    init_grid(1 << minlevel);

    mu = muv;

    run();
}

event init(i = 0) {
    FILE* fp = fopen("NACA0012_processed.gnu", "r");
    points = input_xy(fp);
    distance(d, points);
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);

    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    boundary(all);
}

u.n[left] = dirichlet(U);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[embed] = dirichlet(0.);
// u.t[embed] = neumann(0.);


event properties(i++) {
    foreach_face()
        muv.x[] = fm.x[] * MU;
}

event logfile(i++) {
    fprintf(stderr, "i:%04d t:%.4e dt:%.4e grid:%ld\n", i, t, dt, grid->n);
}

event adapt(i++) {
    adapt_wavelet({cs, u}, (double[]){1e-30, 1e-2, 1e-2}, maxlevel, minlevel);
}

event output(t += 1.0) {
    char name[80];
    static int count = 0;
    sprintf(name, "output/grid_%04d.png", count++);
    view(width = 2000, height = 2000);
    clear();
    squares("u.x", min=0., max=1.5, map=viridis, linear=false);
    save(name);

    FILE* f_ue = fopen("data/ue.dat", "w");
    fprintf(f_ue, "# x y ue\n");
    foreach() {
        if (cs[] > 0. && cs[] < 1.) {
            coord n, b;
            embed_geometry(point, &b, &n);
            fprintf(f_ue, "%.6e %.6e\n", x, embed_interpolate(point, u.x, b));
        }
    }
    fflush(f_ue);
    fclose(f_ue);

    system("gnuplot plot_ue.gp");
}

event end(t = 100.0);