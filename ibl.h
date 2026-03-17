#define is_cutcell (cs[] > 0. && cs[] < 1.)
#define is_fluid (cs[] >= 1.)
#define is_solid (cs[] <= 0.)

Cache cutcells = {0};

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


const double ibl_dt = 1e-5;
scalar ibl_theta[];

#include "tangential_gradients.h"

event init(i = 0) {
    ibl_theta.refine = ibl_theta.prolongation = refine_embed_linear;
    ibl_theta.restriction = restriction_volume_average;
    ibl_theta.dirty = true;

    foreach() {
        if (is_cutcell)
            ibl_theta[] = 1e-10;
        else
            ibl_theta[] = 0.;
    }
}

#if dimension == 2
event update(i++) {
    // H* (Eq. 10)
    double H_star = (Hk < 4.) ?
        1.515 + 0.076 * sq(4. - Hk) / Hk :
        1.515 + 0.040 * sq(Hk - 4.) / Hk;

    // Cf (Eq. 11)
    double Re_Cf2 = (Hk < 7.4) ?
        -0.067 + 0.01977 * sq(7.4 - Hk) / (Hk - 1.) :
        -0.067 + 0.022 * sq(1. - 1.4 / (Hk - 6.));

    double Cf = 2. * Re_Cf2 / Re_theta;

    // Cd (Eq. 12)
    double Re_Cd = (Hk < 4.) ?
        0.207 + 0.00205 * pow(4. - Hk, 5.5) :
        0.207 - 0.003 * sq(Hk - 4.) / (1. + 0.02 * sq(Hk - 4.));

    double Cd = H_star * Re_Cd / (2. * Re_theta);
}
#endif