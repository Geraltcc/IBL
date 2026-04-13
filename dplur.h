#include "run.h"

scalar rho[], E[];
vector w[];

scalar drho[], dE[];
vector dw[];
scalar diag_wall[];

double gammao = 1.4;
int dplur_nrelax = 3;
double CFL_impl = 10.;
double CFL_expl = 0.3;

face vector spec_rad[];

#define SEPS 1e-30

#if EMBED
rho[embed] = neumann(0.);
E[embed] = neumann(0.);
w.n[embed] = dirichlet(0.);
w.t[embed] = neumann(0.);
#endif

void compute_spectral_radius() {
    foreach_face() {
        if (fm.x[] < SEPS) {
            spec_rad.x[] = 0.;
            continue;
        }

        double rhoL = rho[-1], rhoR = rho[0];
        if (rhoL < SEPS && rhoR < SEPS) {
            spec_rad.x[] = 0.;
            continue;
        }
        if (rhoL < SEPS) rhoL = rhoR;
        if (rhoR < SEPS) rhoR = rhoL;

        double uL = w.x[-1]/rhoL, uR = w.x[]/rhoR;
        double w2L = 0., w2R = 0.;
        foreach_dimension() {
            w2L += sq(w.x[-1]/rhoL);
            w2R += sq(w.x[]/rhoR);
        }

        double pL = (gammao - 1.) * (E[-1] - 0.5*rhoL*w2L);
        double pR = (gammao - 1.) * (E[] - 0.5*rhoR*w2R);
        pL = fmax(pL, SEPS); pR = fmax(pR, SEPS);

        double cL = sqrt(gammao * pL / rhoL);
        double cR = sqrt(gammao * pR / rhoR);
        spec_rad.x[] = fmax(fabs(uL) + cL, fabs(uR) + cR);
    }
    boundary({spec_rad});
}

double Roe_flux(const double * right, const double * left, double * f, int len) {
    double rhoL = left[0], rhoR = right[0];
    double EL   = left[1], ER   = right[1];
    double wnL  = left[2], wnR  = right[2];

    if (rhoL < SEPS && rhoR < SEPS) {
        for (int i = 0; i < len; i++) f[i] = 0.;
        return 0.;
    }
    if (rhoL < SEPS) {
        rhoL = rhoR; EL = ER;
        for (int i = 0; i < len; i++) ((double*)left)[i] = right[i];
    }
    if (rhoR < SEPS) {
        rhoR = rhoL; ER = EL;
        for (int i = 0; i < len; i++) ((double*)right)[i] = left[i];
    }

    // --- Primitive variables ---
    double uL = wnL / rhoL, uR = wnR / rhoR;
    double w2L = 0., w2R = 0.;
    for (int i = 2; i < 2 + dimension; i++) {
        w2L += sq(left[i] / rhoL);
        w2R += sq(right[i] / rhoR);
    }
    double pL = (gammao - 1.)*(EL - 0.5*rhoL*w2L);
    double pR = (gammao - 1.)*(ER - 0.5*rhoR*w2R);
    pL = fmax(pL, SEPS);
    pR = fmax(pR, SEPS);
    double HL = (EL + pL)/rhoL;
    double HR = (ER + pR)/rhoR;

    // --- Left/Right physical fluxes ---
    // F = (w_n, u_n(E+p), u_n*w_n + p, u_n*w_t)
    double fL[len], fR[len];
    fL[0] = wnL;
    fL[1] = uL*(EL + pL);
    fL[2] = uL*wnL + pL;
    fR[0] = wnR;
    fR[1] = uR*(ER + pR);
    fR[2] = uR*wnR + pR;
    for (int i = 3; i < 2 + dimension; i++) {
        fL[i] = uL*left[i];
        fR[i] = uR*right[i];
    }

    // --- Roe averages ---
    double sqrL = sqrt(rhoL), sqrR = sqrt(rhoR);
    double denom = sqrL + sqrR;
    double u_roe = (sqrL*uL + sqrR*uR)/denom;
    double H_roe = (sqrL*HL + sqrR*HR)/denom;
    double q2_roe = sq(u_roe);

    // tangential velocities
    double vL = (dimension > 1) ? left[3]/rhoL  : 0.;
    double vR = (dimension > 1) ? right[3]/rhoR : 0.;
    double v_roe = (sqrL*vL + sqrR*vR)/denom;
    q2_roe += sq(v_roe);

    double c2_roe = (gammao - 1.)*(H_roe - 0.5*q2_roe);
    if (c2_roe < SEPS) c2_roe = SEPS;
    double c_roe = sqrt(c2_roe);

    // --- Eigenvalues ---
    double lam[4];
    lam[0] = u_roe - c_roe;  // acoustic-
    lam[1] = u_roe;           // entropy
    lam[2] = u_roe;           // shear
    lam[3] = u_roe + c_roe;  // acoustic+

    // Harten entropy fix
    double eps = 0.1*c_roe;
    for (int k = 0; k < 4; k++) {
        double la = fabs(lam[k]);
        lam[k] = (la >= eps) ? la : (sq(lam[k]) + sq(eps))/(2.*eps);
    }

    // --- Wave strengths ---
    double drho = rhoR - rhoL;
    double du   = uR - uL;
    double dv   = vR - vL;
    double dp   = pR - pL;
    double alpha1 = (dp - sqrt(rhoL*rhoR)*c_roe*du)/(2.*c2_roe);
    double alpha2 = drho - dp/c2_roe;
    double alpha3 = sqrt(rhoL*rhoR)*dv;
    double alpha4 = (dp + sqrt(rhoL*rhoR)*c_roe*du)/(2.*c2_roe);

    // --- Right eigenvectors (ordering: ρ, E, w_n, w_t) ---
    //   r1 = (1,  H̃-ũc̃,  ũ-c̃,  ṽ)
    //   r2 = (1,  ½q̃²,    ũ,     ṽ)
    //   r3 = (0,  ṽ,       0,     1)
    //   r4 = (1,  H̃+ũc̃,  ũ+c̃,  ṽ)

    // --- Dissipation = Σ |λ_k| α_k r_k ---
    double diss[len];
    for (int i = 0; i < len; i++) diss[i] = 0.;
    // Wave 1: acoustic-
    diss[0] += lam[0]*alpha1 * 1.;
    diss[1] += lam[0]*alpha1 * (H_roe - u_roe*c_roe);
    diss[2] += lam[0]*alpha1 * (u_roe - c_roe);
    if (dimension > 1)
        diss[3] += lam[0]*alpha1 * v_roe;
    // Wave 2: entropy
    diss[0] += lam[1]*alpha2 * 1.;
    diss[1] += lam[1]*alpha2 * 0.5*q2_roe;
    diss[2] += lam[1]*alpha2 * u_roe;
    if (dimension > 1)
        diss[3] += lam[1]*alpha2 * v_roe;
    // Wave 3: shear
    if (dimension > 1) {
        diss[1] += lam[2]*alpha3 * v_roe;
        diss[3] += lam[2]*alpha3 * 1.;
    }
    // Wave 4: acoustic+
    diss[0] += lam[3]*alpha4 * 1.;
    diss[1] += lam[3]*alpha4 * (H_roe + u_roe*c_roe);
    diss[2] += lam[3]*alpha4 * (u_roe + c_roe);
    if (dimension > 1)
        diss[3] += lam[3]*alpha4 * v_roe;
    // --- Roe flux = ½(FL + FR) - ½·diss ---
    for (int i = 0; i < len; i++)
        f[i] = 0.5*(fL[i] + fR[i]) - 0.5*diss[i];
    return fmax(fabs(lam[0]), fabs(lam[3]));
}

void mirror_state(const double * fluid, double * ghost, const double * n, int len) {
    ghost[0] = fluid[0]; // rho
    ghost[1] = fluid[1]; // E

    double wdotn = 0.;
    for (int d = 0; d < dimension; d++)
        wdotn += fluid[2 + d] * n[d];
    for (int d = 0; d < dimension; d++)
        ghost[2 + d] = fluid[2 + d] - 2.*wdotn*n[d];
}

event defaults(i = 0) {
    CFL = CFL_expl;

    for (scalar s in {rho, E, w})
        s.gradient = minmod2;

#if TREE
    for (scalar s in {rho, E, w, drho, dE, dw}) {
#if EMBED
        s.refine = s.prolongation = refine_embed_linear;
        s.restriction = restriction_volume_average;
        s.depends = list_add(s.depends, cs);
#else
        s.refine = s.prolongation = refine_linear;
#endif
        s.dirty = true;
    }

    drho.nodump = dE.nodump = true;
    foreach_dimension()
        dw.x.nodump = true;
#endif
}

event init(i = 0) {
#if EMBED
    fractions_cleanup(cs, fs);
#endif
    boundary({rho, E, w});
    dtmax = DT;
    dt = dtnext(1e-6);
}

double dtmax;

event flux_compute(i++, last) {
    // Zero solid cells
    foreach()
        if (cs[] <= 0.)
            for (scalar s in {rho, E, w})
                s[] = 0.;

    // ---- Precompute embed normals ----
#if EMBED
    vector n_embed = new vector;
    foreach() {
        foreach_dimension()
            n_embed.x[] = 0.;
        if (cs[] > 0. && cs[] < 1.) {
            coord n, b;
            embed_geometry(point, &b, &n);
            foreach_dimension()
                n_embed.x[] = n.x;
        }
    }
    boundary((scalar *){n_embed});
#endif

    // ---- MUSCL gradients (minmod limiter) ----
    int len = 2 + dimension;   // ρ, E, w.x, [w.y]
    scalar * conserved = {rho, E, w};
    
    vector * slopes = NULL;
    for (scalar s in conserved) {
        vector slope = new vector;
        foreach_dimension() {
            slope.x.gradient = minmod2;
#if TREE && EMBED
            slope.x.prolongation = refine_embed_linear;
#elif TREE
            slope.x.prolongation = refine_linear;
#endif
        }
        slopes = vectors_append(slopes, slope);
    }

    foreach() {
        scalar s; vector g;
        for (s, g in conserved, slopes) {
            foreach_dimension() {
#if EMBED
                if (!fs.x[] && !fs.x[1])
                    g.x[] = 0.;
                else if (!fs.x[])
                    g.x[] = s.gradient(s[], s[], s[1])/Delta;
                else if (!fs.x[1])
                    g.x[] = s.gradient(s[-1], s[], s[])/Delta;
                else
#endif
                    g.x[] = s.gradient(s[-1], s[], s[1])/Delta;
            }
        }
    }
    boundary((scalar *)slopes);

    // ---- Allocate face fluxes ----
    vector * lflux = NULL;
    for (scalar s in conserved) {
        vector f1 = new face vector;
        lflux = vectors_append(lflux, f1);
    }
    for (vector fl in lflux)
        foreach_face()
            fl.x[] = 0.;

    double dt_ref = HUGE;
    // ---- Compute face fluxes ----
    foreach_face(reduction(min:dt_ref)) {
        if (fm.x[] < SEPS) continue;
        double r[len], l[len], f[len];
        double dx = Delta/2.;
        // MUSCL reconstruction
        int idx = 0;
        scalar s; vector g;
        for (s, g in conserved, slopes) {
            r[idx] = s[]   - dx*g.x[];
            l[idx] = s[-1] + dx*g.x[-1];
            idx++;
        }

        double face_a;
        bool right_ok = (cs[]   > 0.);
        bool left_ok  = (cs[-1] > 0.);
        if (!right_ok && !left_ok) {
            for (int k = 0; k < len; k++) f[k] = 0.;
        }
#if EMBED
        else if (!right_ok && left_ok) {
            double n_wall[dimension], nmag2 = 0.;
            n_wall[0] = n_embed.x[];
#if dimension > 1
            n_wall[1] = n_embed.y[];
#endif
            for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
            if (nmag2 < SEPS) {
                n_wall[0] = n_embed.x[-1];
#if dimension > 1
                n_wall[1] = n_embed.y[-1];
#endif
                nmag2 = 0.;
                for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
            }
            if (nmag2 < SEPS) {
                n_wall[0] = 1.;
                for (int d = 1; d < dimension; d++) n_wall[d] = 0.;
                nmag2 = 1.;
            }
            double nmag = sqrt(nmag2);
            for (int d = 0; d < dimension; d++) n_wall[d] /= nmag;
            mirror_state(l, r, n_wall, len);
            Roe_flux(r, l, f, len);
        }
        else if (right_ok && !left_ok) {
            double n_wall[dimension], nmag2 = 0.;
            n_wall[0] = n_embed.x[-1];
#if dimension > 1
            n_wall[1] = n_embed.y[-1];
#endif
            for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
            if (nmag2 < SEPS) {
                n_wall[0] = n_embed.x[];
#if dimension > 1
                n_wall[1] = n_embed.y[];
#endif
                nmag2 = 0.;
                for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
            }
            if (nmag2 < SEPS) {
                n_wall[0] = 1.;
                for (int d = 1; d < dimension; d++) n_wall[d] = 0.;
                nmag2 = 1.;
            }
            double nmag = sqrt(nmag2);
            for (int d = 0; d < dimension; d++) n_wall[d] /= nmag;
            mirror_state(r, l, n_wall, len);
            Roe_flux(r, l, f, len);
        }
#endif // EMBED
        else {
            face_a = Roe_flux(r, l, f, len);
        }

        if (right_ok && left_ok && face_a > 0.) {
            double dtc = Delta / face_a;
            if (dtc < dt_ref) dt_ref = dtc;
        }

        idx = 0;
        for (vector fl in lflux)
            fl.x[] = fm.x[] * f[idx++];
    }

    // ---- Divergence → residual stored in (drho, dE, dw) ----
    scalar * res = {drho, dE, dw};
    foreach() {
        if (cs[] <= 0.) {
            for (scalar ds in res)
                ds[] = 0.;
            continue;
        }
        scalar ds; vector fl;
        for (ds, fl in res, lflux) {
            ds[] = 0.;
            foreach_dimension()
                ds[] += (fl.x[] - fl.x[1])/(cm[]*Delta);
        }
    }

    // ---- Embedded surface pressure flux ----
#if EMBED
    foreach()
        diag_wall[] = 0.;
    foreach() {
        if (cs[] <= 0. || cs[] >= 1.) continue;
        coord n, b;
        double area = embed_geometry(point, &b, &n);
        double rho_b = embed_interpolate(point, rho, b);
        double E_b   = embed_interpolate(point, E, b);
        double w2_b  = 0.;
        foreach_dimension() {
            double wb = embed_interpolate(point, w.x, b);
            w2_b += sq(wb);
        }
        double p_wall = (gammao - 1.)*(E_b - 0.5*w2_b/fmax(rho_b, SEPS));
        p_wall = fmax(p_wall, 0.);
        double cm_safe = fmax(cs[], 1e-2);
        foreach_dimension()
            dw.x[] -= p_wall*area*n.x / (cm_safe*Delta);

        // Wall face spectral radius for DPLUR diagonal
        double c_wall = sqrt(gammao*p_wall/fmax(rho_b, SEPS));
        diag_wall[] = 0.5*area*(c_wall) / (cm_safe*Delta);
    }
    boundary({diag_wall});
#endif

    // ---- Cleanup ----
#if EMBED
    delete((scalar *){n_embed});
    #endif
    delete((scalar *)slopes);  free(slopes);
    delete((scalar *)lflux);   free(lflux);
    boundary(res);

    // dt = dtnext(fmin(CFL_impl * dt_ref, dtmax));
    dt = dtnext(1e-6);
}

// ============================================================
//  Event: dplur_step — DPLUR relaxation + state update
// ============================================================

event dplur_step(i++, last)
{
    compute_spectral_radius();
    scalar * delta = {drho, dE, dw};
    // Save residual
    scalar * res_save = list_clone(delta);
    foreach() {
        scalar ds, rs;
        for (ds, rs in delta, res_save)
            rs[] = ds[];
    }
    // Iteration 0: δU = D⁻¹ · R
    foreach() {
        if (cs[] <= 0.) {
            for (scalar ds in delta) ds[] = 0.;
            continue;
        }
        double inv_cmd = 1./(cm[]*Delta);
        double diag = 1./dt;
        foreach_dimension()
            diag += 0.5*(fm.x[]*spec_rad.x[]
                       + fm.x[1]*spec_rad.x[1]) * inv_cmd;

#if EMBED
        diag += diag_wall[];
#endif
        double inv_d = 1./diag;
        for (scalar ds in delta)
            ds[] *= inv_d;
    }
    // Iterations 1 .. nrelax-1
    for (int k = 1; k < dplur_nrelax; k++) {
        boundary(delta);
        foreach() {
            if (cs[] <= 0.) continue;
            double inv_cmd = 1./(cm[]*Delta);
            double diag = 1./dt;
            foreach_dimension()
                diag += 0.5*(fm.x[]*spec_rad.x[]
                           + fm.x[1]*spec_rad.x[1]) * inv_cmd;

                diag += diag_wall[];
            double inv_d = 1./diag;
            scalar ds, rs;
            for (ds, rs in delta, res_save) {
                double nbr = 0.;
                foreach_dimension()
                    nbr += 0.5*(fm.x[]*spec_rad.x[]*ds[-1]
                              + fm.x[1]*spec_rad.x[1]*ds[1]) * inv_cmd;
                ds[] = (rs[] + nbr)*inv_d;
            }
        }
    }
    // Apply: U^{n+1} = U^n + δU
    foreach() {
        if (cs[] <= 0.) continue;
        rho[] += drho[];
        E[]   += dE[];
        foreach_dimension()
            w.x[] += dw.x[];
        // Positivity
        if (rho[] < SEPS) rho[] = SEPS;
        double w2 = 0.;
        foreach_dimension()
            w2 += sq(w.x[]);
        double pc = (gammao - 1.)*(E[] - 0.5*w2/rho[]);
        if (pc < SEPS)
            E[] = SEPS/(gammao - 1.) + 0.5*w2/rho[];
    }
    delete(res_save);
    free(res_save);
    boundary({rho, E, w});
}

#if TREE
event adapt(i++, last)
{
#if EMBED
    fractions_cleanup(cs, fs);
#endif
}
#endif
