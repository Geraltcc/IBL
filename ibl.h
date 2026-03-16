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

void cutcell_tangential_gradient(scalar f, vector g) {
    foreach() 
        foreach_dimension()
            g.x[] = 0.;
    
    foreach_cache(cutcells) {
        coord n, b, grad = {0.};
        embed_geometry(point, &b, &n);

        foreach_dimension() {
            bool has_p = (cs[1] > 0. && cs[1] < 1.);
            bool has_m = (cs[-1] > 0. && cs[-1] < 1.);

            if (has_p && has_m) 
                grad.x = (f[1] - f[-1]) / (2.*Delta);
            else if (has_p)
                grad.x = (f[1] - f[]) / Delta;
            else if (has_m)
                grad.x = (f[] - f[-1]) / Delta;
        }

        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n.x;

        foreach_dimension()
            g.x[] = grad.x - gn * n.x;
    }
}

void cutcell_tangential_gradient_extend(scalar f, vector g) {
    scalar f_ext[];
    foreach() {
        f_ext[] = f[];
        foreach_dimension()
            g.x[] = 0.;

        if (cs[] > 0.)
            continue;

        double sum = 0.;
        int count = 0.;
        foreach_dimension() {
            if (cs[1] > 0. && cs[1] < 1.) { sum += f[1]; count++; }
            if (cs[-1] > 0. && cs[-1] < 1.) { sum += f[-1]; count++; }
        }
        if (count > 0)
            f_ext[] = sum / count;
    }

    foreach_cache(cutcells) {
        coord n, b;
        embed_geometry(point, &b, &n);

        coord grad = {0.};
        foreach_dimension()
            grad.x = center_gradient(f_ext);

        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n.x;

        foreach_dimension()
            g.x[] = grad.x - gn * n.x;
    }
    delete({f_ext});
}

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

/*
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
*/