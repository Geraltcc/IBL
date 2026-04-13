#include "grid/octree.h"
#include "embed.h"
#include "distance.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "cmap.h"

#define maxlevel 8
#define minlevel 6

face vector muv[];

int main() {
    L0 = 100.0;
    size(L0);
    origin(-L0/2., -L0/2., -L0/2.);

    init_grid(1 << minlevel);

    mu = muv;

    run();
}

scalar d[];
coord* points = NULL;

event init(i = 0) {
    FILE* fp = fopen("airfoil/CHN-F1-binary.stl", "r");
    points = input_stl (fp);

    distance(d, points);
    while (adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);

    foreach()
        d[] *= -1.;

    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    foreach_face() 
        muv.x[] = fm.x[] * 1e-3;

    foreach() {
        u.x[] = cs[];
    }

    view();
    clear();
    squares("cs", map=viridis, linear=false);
    save("output/cs.png");
}

u.n[left] = dirichlet(1.);
p[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

event logfile(i++) {
    fprintf(stderr, "i: %04d t: %.4e dt: %.4e grid: %ld\n", i, t, dt, grid->n);
}

event adapt(i++) {
    adapt_wavelet({cs, u}, (double[]){1e-30, 1e-2}, maxlevel, minlevel);
}

event properties(i++) {
    foreach_face() 
        muv.x[] = fm.x[] * 1e-3;
}

event output(i += 100) {
    char name[80];
    sprintf(name, "output/grid_%04d.png", i);
    view(width = 2000, height = 2000);
    clear();
    squares("u.x", min=0., max=1.5, map=viridis, linear=false);
    save(name);
}