#include "grid/quadtree.h"
#include "view.h"
#include "cmap.h"

#define x_offset -1.0

int main() {
    restore("restart");
    boundary(all);
    view(width = 2000, height = 2000);
    clear();
    squares("w.x", min=0., max=1.5, map=viridis, linear=false);
    save("output/w.x.png");
}
