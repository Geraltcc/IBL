/**
# Implicit compressible Euler solver (Point-Implicit + Jacobi)

Follows the predictor-corrector.h / conservation.h structural pattern:
- run() provides the main time loop with dt managed via dtnext()
- update_conservation() computes Roe residual and returns dtmax
- advance_implicit() applies point-implicit + Jacobi relaxation

Usage: identical to mycompressible_nofr.h.
*/

#include "utils.h"

#define SEPS 1e-30
#define CS_FLOOR 0.1

// ==================== Conserved fields ====================
scalar rho[], E[];
vector w[];
scalar * scalars = {rho, E};
vector * myvectors = {w};
double gammao = 1.4;

// ==================== Solver parameters ====================
double CFL_impl = 50.;  // implicit CFL (50–500, ramp during convergence)
int    nrelax   = 3;     // Jacobi relaxation iterations

// ==================== Internal fields ====================
scalar * evolving;
double dt = 0.;
face vector spec_rad[];
scalar diag_impl[];

// ==================== Embedded BCs ====================
#if EMBED
rho[embed]  = neumann(0.);
E[embed]    = neumann(0.);
w.n[embed]  = dirichlet(0.);
w.t[embed]  = neumann(0.);
#endif

// ==================== Spectral radius ====================
static void compute_spectral_radius()
{
  foreach_face() {
    if (fm.x[] < SEPS) { spec_rad.x[] = 0.; continue; }
    double rhoL = rho[-1], rhoR = rho[];
    if (rhoL < SEPS && rhoR < SEPS) { spec_rad.x[] = 0.; continue; }
    if (rhoL < SEPS) rhoL = rhoR;
    if (rhoR < SEPS) rhoR = rhoL;
    double uL = w.x[-1]/rhoL, uR = w.x[]/rhoR;
    double w2L = 0., w2R = 0.;
    foreach_dimension() {
      w2L += sq(w.x[-1]/rhoL);
      w2R += sq(w.x[]/rhoR);
    }
    double pL = fmax((gammao-1.)*(E[-1] - 0.5*rhoL*w2L), SEPS);
    double pR = fmax((gammao-1.)*(E[]   - 0.5*rhoR*w2R), SEPS);
    spec_rad.x[] = fmax(fabs(uL) + sqrt(gammao*pL/rhoL),
                        fabs(uR) + sqrt(gammao*pR/rhoR));
  }
  boundary((scalar *){spec_rad});
}

// ==================== Roe flux ====================
static void mirror_state (const double * fluid, double * ghost,
                          const double * n, int len)
{
  ghost[0] = fluid[0];
  ghost[1] = fluid[1];
  double wdotn = 0.;
  for (int d = 0; d < dimension; d++)
    wdotn += fluid[2+d]*n[d];
  for (int d = 0; d < dimension; d++)
    ghost[2+d] = fluid[2+d] - 2.*wdotn*n[d];
}

static double Roe_flux (const double * right, const double * left,
                        double * f, int len)
{
  double rhoL = left[0],  rhoR = right[0];
  double EL   = left[1],  ER   = right[1];
  double wnL  = left[2],  wnR  = right[2];
  if (rhoL < SEPS && rhoR < SEPS) {
    for (int i = 0; i < len; i++) f[i] = 0.;
    return 0.;
  }
  if (rhoL < SEPS) { rhoL = rhoR; EL = ER;
    for (int i = 0; i < len; i++) ((double*)left)[i] = right[i]; }
  if (rhoR < SEPS) { rhoR = rhoL; ER = EL;
    for (int i = 0; i < len; i++) ((double*)right)[i] = left[i]; }

  double uL = wnL/rhoL, uR = wnR/rhoR;
  double w2L = 0., w2R = 0.;
  for (int i = 2; i < 2+dimension; i++) {
    w2L += sq(left[i]/rhoL);
    w2R += sq(right[i]/rhoR);
  }
  double pL = fmax((gammao-1.)*(EL - 0.5*rhoL*w2L), SEPS);
  double pR = fmax((gammao-1.)*(ER - 0.5*rhoR*w2R), SEPS);
  double HL = (EL+pL)/rhoL, HR = (ER+pR)/rhoR;

  double fL[len], fR[len];
  fL[0] = wnL;           fR[0] = wnR;
  fL[1] = uL*(EL+pL);    fR[1] = uR*(ER+pR);
  fL[2] = uL*wnL + pL;   fR[2] = uR*wnR + pR;
  for (int i = 3; i < 2+dimension; i++) {
    fL[i] = uL*left[i];  fR[i] = uR*right[i];
  }

  double sqrL = sqrt(rhoL), sqrR = sqrt(rhoR), den = sqrL+sqrR;
  double u_roe = (sqrL*uL + sqrR*uR)/den;
  double H_roe = (sqrL*HL + sqrR*HR)/den;
  double q2 = sq(u_roe);
  double vL = (dimension > 1) ? left[3]/rhoL  : 0.;
  double vR = (dimension > 1) ? right[3]/rhoR : 0.;
  double v_roe = (sqrL*vL + sqrR*vR)/den;
  q2 += sq(v_roe);

  double c2 = fmax((gammao-1.)*(H_roe - 0.5*q2), SEPS);
  double c  = sqrt(c2);

  double lam[4] = { u_roe-c, u_roe, u_roe, u_roe+c };
  double eps_h = 0.1*c;
  for (int k = 0; k < 4; k++) {
    double la = fabs(lam[k]);
    lam[k] = (la >= eps_h) ? la : (sq(lam[k])+sq(eps_h))/(2.*eps_h);
  }

  double drho = rhoR-rhoL, du = uR-uL, dv = vR-vL, dp = pR-pL;
  double rhoG = sqrt(rhoL*rhoR);
  double a1 = (dp - rhoG*c*du)/(2.*c2);
  double a2 = drho - dp/c2;
  double a3 = rhoG*dv;
  double a4 = (dp + rhoG*c*du)/(2.*c2);

  double diss[len];
  for (int i = 0; i < len; i++) diss[i] = 0.;
  diss[0] += lam[0]*a1;
  diss[1] += lam[0]*a1*(H_roe - u_roe*c);
  diss[2] += lam[0]*a1*(u_roe - c);
  if (dimension > 1) diss[3] += lam[0]*a1*v_roe;

  diss[0] += lam[1]*a2;
  diss[1] += lam[1]*a2*0.5*q2;
  diss[2] += lam[1]*a2*u_roe;
  if (dimension > 1) diss[3] += lam[1]*a2*v_roe;

  if (dimension > 1) { diss[1] += lam[2]*a3*v_roe; diss[3] += lam[2]*a3; }

  diss[0] += lam[3]*a4;
  diss[1] += lam[3]*a4*(H_roe + u_roe*c);
  diss[2] += lam[3]*a4*(u_roe + c);
  if (dimension > 1) diss[3] += lam[3]*a4*v_roe;

  for (int i = 0; i < len; i++)
    f[i] = 0.5*(fL[i] + fR[i]) - 0.5*diss[i];
  return fmax(fabs(lam[0]), fabs(lam[3]));
}

// ==================== Flux computation & residual ====================
static double update_conservation (scalar * conserved, scalar * updates,
                                   double dtmax)
{
  foreach()
    if (cs[] <= 0.)
      for (scalar s in conserved)
        s[] = 0.;

#if EMBED
  vector n_embed = new vector;
  foreach() {
    foreach_dimension() n_embed.x[] = 0.;
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      embed_geometry(point, &b, &n);
      foreach_dimension() n_embed.x[] = n.x;
    }
  }
  boundary((scalar *){n_embed});
#endif

  int len = list_len(conserved);
  int scalars_len = list_len(scalars);

  // MUSCL gradients
  vector * slopes = NULL;
  for (scalar s in conserved) {
    vector slope = new vector;
    foreach_dimension() {
      slope.x.gradient = zero;
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
    for (s,g in conserved,slopes) {
      foreach_dimension() {
#if EMBED
        if (!fs.x[] && !fs.x[1])       g.x[] = 0.;
        else if (!fs.x[])              g.x[] = minmod2(s[], s[], s[1])/Delta;
        else if (!fs.x[1])             g.x[] = minmod2(s[-1], s[], s[])/Delta;
        else
#endif
          g.x[] = minmod2(s[-1], s[], s[1])/Delta;
      }
    }
  }
  boundary((scalar *)slopes);

  // Split conserved list into scalars/vectors for correct face rotation
  scalar * sc_list   = list_copy(conserved);
  if (sc_list) sc_list[scalars_len].i = -1;
  vector * vec_list  = vectors_from_scalars(&conserved[scalars_len]);

  vector * sc_slopes = vectors_copy(slopes);
  if (sc_slopes) sc_slopes[scalars_len] = (vector){{-1}};
  tensor * vec_slopes = tensors_from_vectors(&slopes[scalars_len]);

  // Face fluxes
  vector * lflux = NULL;
  for (scalar s in conserved) {
    vector f1 = new face vector;
    lflux = vectors_append(lflux, f1);
  }

  vector * sc_fluxes = vectors_copy(lflux);
  if (sc_fluxes) sc_fluxes[scalars_len] = (vector){{-1}};
  tensor * vec_fluxes = tensors_from_vectors(&lflux[scalars_len]);

  for (vector fl in lflux)
    foreach_face() fl.x[] = 0.;

  double dt_ref = HUGE;
  foreach_face(reduction(min:dt_ref)) {
    if (fm.x[] < SEPS) continue;
    double r[len], l[len], f[len];
    double dx = Delta/2.;
    int i = 0;
    scalar s; vector g;
    for (s,g in sc_list, sc_slopes) {
      r[i]   = s[]   - dx*g.x[];
      l[i++] = s[-1] + dx*g.x[-1];
    }
    vector v; tensor t;
    for (v,t in vec_list, vec_slopes) {
      r[i]   = v.x[]   - dx*t.x.x[];
      l[i++] = v.x[-1] + dx*t.x.x[-1];
#if dimension > 1
      r[i]   = v.y[]   - dx*t.y.x[];
      l[i++] = v.y[-1] + dx*t.y.x[-1];
#endif
    }

    bool right_ok = (cs[] > 0.), left_ok = (cs[-1] > 0.);
    double face_a = 0.;

    if (!right_ok && !left_ok) {
      for (int k = 0; k < len; k++) f[k] = 0.;
    }
#if EMBED
    else if (left_ok && !right_ok) {
      double n_wall[dimension];
      n_wall[0] = n_embed.x[];
#if dimension > 1
      n_wall[1] = n_embed.y[];
#endif
      double nmag2 = 0.;
      for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
      if (nmag2 < SEPS) {
        n_wall[0] = n_embed.x[-1];
#if dimension > 1
        n_wall[1] = n_embed.y[-1];
#endif
        nmag2 = 0.;
        for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
      }
      if (nmag2 < SEPS) { n_wall[0] = 1.; nmag2 = 1.; }
      double nmag = sqrt(nmag2);
      for (int d = 0; d < dimension; d++) n_wall[d] /= nmag;
      mirror_state(l, r, n_wall, len);
      Roe_flux(r, l, f, len);
    }
    else if (right_ok && !left_ok) {
      double n_wall[dimension];
      n_wall[0] = n_embed.x[-1];
#if dimension > 1
      n_wall[1] = n_embed.y[-1];
#endif
      double nmag2 = 0.;
      for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
      if (nmag2 < SEPS) {
        n_wall[0] = n_embed.x[];
#if dimension > 1
        n_wall[1] = n_embed.y[];
#endif
        nmag2 = 0.;
        for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);
      }
      if (nmag2 < SEPS) { n_wall[0] = 1.; nmag2 = 1.; }
      double nmag = sqrt(nmag2);
      for (int d = 0; d < dimension; d++) n_wall[d] /= nmag;
      mirror_state(r, l, n_wall, len);
      Roe_flux(r, l, f, len);
    }
#endif
    else {
      face_a = Roe_flux(r, l, f, len);
      if (face_a > 0. && cs[] >= CS_FLOOR && cs[-1] >= CS_FLOOR) {
        double dtc = Delta / face_a;
        if (dtc < dt_ref) dt_ref = dtc;
      }
    }

    i = 0;
    for (vector fs in sc_fluxes)  fs.x[] = fm.x[]*f[i++];
    for (tensor fv in vec_fluxes) {
      fv.x.x[] = fm.x[]*f[i++];
#if dimension > 1
      fv.y.x[] = fm.x[]*f[i++];
#endif
    }
  }

  // Divergence → updates
  foreach() {
    if (cs[] <= 0.) {
      for (scalar ds in updates) ds[] = 0.;
      continue;
    }
    scalar ds; vector fl;
    for (ds,fl in updates,lflux) {
      ds[] = 0.;
      foreach_dimension()
        ds[] += (fl.x[] - fl.x[1])/(cm[]*Delta);
    }
  }

  // Embedded wall pressure
#if EMBED
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
    double p_wall = fmax((gammao-1.)*(E_b - 0.5*w2_b/fmax(rho_b,SEPS)), 0.);
    double cm_safe = fmax(cs[], 1e-2);
    int idx = 0;
    for (scalar ds in updates) {
      if (idx >= scalars_len && idx < scalars_len + dimension) {
        int d = idx - scalars_len;
        double nc = (d==0) ? n.x :
#if dimension > 1
                    (d==1) ? n.y :
#endif
                    0.;
        ds[] -= p_wall*area*nc / (cm_safe*Delta);
      }
      idx++;
    }
  }
#endif

  // Spectral radius + implicit diagonal
  compute_spectral_radius();
  foreach() {
    diag_impl[] = 0.;
    if (cs[] <= 0.) continue;
    double inv_cmd = 1./(cm[]*Delta);
    foreach_dimension()
      diag_impl[] += 0.5*(fm.x[]*spec_rad.x[]
                        + fm.x[1]*spec_rad.x[1]) * inv_cmd;
#if EMBED
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry(point, &b, &n);
      double rho_b = embed_interpolate(point, rho, b);
      double E_b   = embed_interpolate(point, E, b);
      double w2_b  = 0.;
      foreach_dimension() {
        double wb = embed_interpolate(point, w.x, b);
        w2_b += sq(wb);
      }
      double p_w = fmax((gammao-1.)*(E_b - 0.5*w2_b/fmax(rho_b,SEPS)), SEPS);
      double c_w = sqrt(gammao*p_w/fmax(rho_b,SEPS));
      diag_impl[] += 0.5*area*c_w / (fmax(cs[],1e-2)*Delta);
    }
#endif
  }
  boundary({diag_impl});

  // Cleanup
#if EMBED
  delete((scalar *){n_embed});
#endif
  free(sc_list); free(vec_list);
  free(sc_slopes); free(vec_slopes);
  free(sc_fluxes); free(vec_fluxes);
  delete((scalar *)slopes); free(slopes);
  delete((scalar *)lflux);  free(lflux);

  return fmin(CFL_impl * dt_ref, dtmax);
}

// ==================== Point-Implicit + Jacobi advance ====================
static void advance_implicit (scalar * output, scalar * input,
                              scalar * updates, double dt_step)
{
  scalar * res_save = NULL;
  if (nrelax > 1) {
    res_save = list_clone(updates);
    foreach() {
      scalar u, r;
      for (u,r in updates, res_save)
        r[] = u[];
    }
  }

  // Iteration 0: ΔU = R / (1/dt + D)
  foreach() {
    if (cs[] <= 0.) {
      for (scalar u in updates) u[] = 0.;
      continue;
    }
    double diag = 1./dt_step + diag_impl[];
    for (scalar u in updates)
      u[] /= diag;
  }

  // Jacobi iterations 1 .. nrelax-1
  for (int k = 1; k < nrelax; k++) {
    boundary(updates);
    foreach() {
      if (cs[] <= 0.) continue;
      double inv_cmd = 1./(cm[]*Delta);
      double diag = 1./dt_step + diag_impl[];
      scalar u, r;
      for (u,r in updates, res_save) {
        double nbr = 0.;
        foreach_dimension()
          nbr += 0.5*(fm.x[]*spec_rad.x[]*u[-1]
                    + fm.x[1]*spec_rad.x[1]*u[1]) * inv_cmd;
        u[] = (r[] + nbr) / diag;
      }
    }
  }

  if (res_save) { delete(res_save); free(res_save); }

  // Apply: output = input + ΔU, with positivity enforcement
  foreach() {
    if (cs[] <= 0.) continue;
    scalar o, i, u;
    for (o,i,u in output, input, updates)
      o[] = i[] + u[];
    rho[] = fmax(rho[], SEPS);
    double w2 = 0.;
    foreach_dimension() w2 += sq(w.x[]);
    double pc = (gammao-1.)*(E[] - 0.5*w2/rho[]);
    if (pc < SEPS)
      E[] = SEPS/(gammao-1.) + 0.5*w2/rho[];
  }
  boundary(output);
}

// ==================== Events ====================
event defaults (i = 0)
{
  evolving = list_concat(scalars, (scalar *)myvectors);

  for (scalar s in evolving)
    s.gradient = minmod2;

  if (CFL > 1./dimension)
    CFL = 1./dimension;

#if TREE
  for (scalar s in evolving) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
    s.restriction = restriction_volume_average;
    s.depends = list_add(s.depends, cs);
#else
    s.refine = s.prolongation = refine_linear;
#endif
    s.dirty = true;
  }
#endif
}

event cleanup (i = end) free(evolving);
event init (i = 0);

// ==================== Main time loop ====================
trace
void run()
{
  t = 0., iter = 0;
  dimensional (dt == DT);
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {
    scalar * updates = list_clone(evolving);
    double dtmax = update_conservation(evolving, updates, DT);
    dt = dtnext(fmax(dtmax, 1e-6));

    advance_implicit(evolving, evolving, updates, dt);

    delete(updates);
    free(updates);
    update_perf();
    iter = inext, t = tnext;
  }
  timer_print(perf.gt, iter, perf.tnc);
  free_grid();
}