#include "grid/quadtree.h"
#include "embed.h"
#include "ibl.h"
#include "view.h"
#include "cmap.h"

scalar f[];

int main() {
    L0 = 1.;
    size(L0);

    origin(-L0/2., -L0/2.);

    init_grid(1 << 6);

    solid(cs, fs, sq(x) + sq(y) - sq(0.25));
    build_cutcell_cache();

    foreach() 
        f[] = 0.;
    
    double sum = 1.;
    foreach_cache(cutcells) {
        f[] += sum;
        sum += 1.;
    }


    view();
    clear();
    squares("f", min=0., max=sum, map=viridis, linear=false);
    save("output/f.png");
}
