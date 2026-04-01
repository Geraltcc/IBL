#define is_cutcell (cs[] > 0. && cs[] < 1.)
#define is_fluid (cs[] >= 1.)
#define is_solid (cs[] <= 0.)

Cache cutcells = {0};

bool is_nested(Point point) {
    bool flag = true;
    foreach_neighbor(1) {
        if (is_solid) {
            flag = false;
            break;
        }
    }
    return flag;
}

void build_cutcell_cache() {
    free(cutcells.p);
    cutcells.p = NULL;
    cutcells.n = cutcells.nm = 0;

    // Pass 1: parallel count
    int count = 0;
    foreach(reduction(+:count))
        if (is_cutcell)
            count++;
    
    // allocate memory
    cutcells.nm = count;
    cutcells.p = (Index *)malloc(count * sizeof(Index));
    cutcells.n = 0;

    // Pass 2: parallel fill
    foreach() {
        if (is_cutcell) {
            int idx = __atomic_fetch_add(&cutcells.n, 1, __ATOMIC_RELAXED);
            cutcells.p[idx].i = point.i;
    #if dimension > 1
            cutcells.p[idx].j = point.j;
    #endif
    #if dimension > 2
            cutcells.p[idx].k = point.k;
    #endif
            cutcells.p[idx].level = point.level;
            cutcells.p[idx].flags = 0;
        }
    }
}

#include "tangential_gradients.h"
