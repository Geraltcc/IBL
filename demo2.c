// #include "grid/quadtree.h"
#include "grid/octree.h"
// #include "embed.h"
// #include "ibl.h"
#include "view.h"
#include "cmap.h"

scalar f[];

int main() {
    origin(-L0/2., -L0/2.);
    init_grid(64);

    // solid (cs, fs, sq(x) + sq(y) - sq(0.25));

    // build_cutcell_cache();

    foreach()
        f[] = 0.;

    int count = 0;
    foreach_point(0., 0., reduction(+:count)) {
        foreach_neighbor(2) {
            count++;
        }
    }
    fprintf(stderr, "count: %d\n", count);

    view();
    clear();
    squares("f", min = 0., max = 1., map = viridis, linear = false);
    save("output/f.png");
}