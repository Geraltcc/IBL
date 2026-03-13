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
        if (cs[] > 0. && cs[] < 1.)
            count++;
    
    // allocate memory
    cutcells.nm = count;
    cutcells.p = (Index *)malloc(count * sizeof(Index));
    cutcells.n = 0;

    // Pass 2: parallel fill
    foreach() {
        if (cs[] > 0. && cs[] < 1.) {
            int idx = __atomic_fetch_add(&cutcells.n, 1, __ATOMIC_RELAXED);
            cutcells.p[idx].i = point.i;
    #if dimension >= 2
            cutcells.p[idx].j = point.j;
    #endif
    #if dimension >= 3
            cutcells.p[idx].k = point.k;
    #endif
            cutcells.p[idx].level = point.level;
            cutcells.p[idx].flags = 0;
        }
    }
}

const double ibl_dt = 1e-5;
scalar ibl_theta[];

void cutcell_tangential_gradient(scalar f, vector g) {

    foreach_cache(cutcells) {
        coord n, b, grad = {0.};
        embed_geometry(point, &b, &n);

        foreach_dimension() {
            double fp = (cs[1] > 0.) ? ((cs[1] < 1.) ? f[1] : f[]) : 0.;
            double fm = (cs[-1] > 0.) ? ((cs[-1] < 1.) ? f[-1] : f[]) : 0.;

            if (fs.x[] && fs.x[1])
                grad.x = (fp - fm) / (2. * Delta);
            else if (fs.x[1])
                grad.x = (fp - f[]) / Delta;
            else if (fs.x[])
                grad.x = (f[] - fm) / Delta;
        }


        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n.x;

        foreach_dimension()
            g.x[] = grad.x - gn * n.x;

        foreach_dimension()
            g.x[] = grad.x;
    }
}

event init(i = 0) {
    ibl_theta.refine = ibl_theta.prolongation = refine_embed_linear;
    ibl_theta.restriction = restriction_volume_average;
    ibl_theta.dirty = true;

    foreach() {
        if (cs[] > 0. && cs[] < 1.)
            ibl_theta[] = 1e-10;
        else
            ibl_theta[] = 0.;
    }
}

event update(i++) {
    
}
