/**
# A generic solver for systems of conservation laws

Using the ideas of [Kurganov and Tadmor,
2000](references.bib#kurganov2000) it is possible to write a generic
solver for systems of conservation laws of the form
$$
\partial_t\left(\begin{array}{c}
    s_i\\
    \mathbf{v}_j\\
 \end{array}\right) + \nabla\cdot\left(\begin{array}{c}
    \mathbf{F}_i\\
    \mathbf{T}_j\\
 \end{array}\right) = 0
$$
where $s_i$ is a list of scalar fields, $\mathbf{v}_j$ a list of
vector fields and $\mathbf{F}_i$, $\mathbf{T}_j$ are the corresponding
vector (resp. tensor) fluxes. 

Note that the [Saint-Venant solver](saint-venant.h) is a particular
case of this generic algorithm.

The user must provide the lists of conserved scalar and vector fields
*/
double my_embed_interpolate (Point point, scalar s, coord b)
{
  int i = sign(b.x), j = sign(b.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j])
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(b.x)) + s[i]*fabs(b.x))*(1. - fabs(b.y)) + 
	    (s[0,j]*(1. - fabs(b.x)) + s[i,j]*fabs(b.x))*fabs(b.y));
#else // dimension == 3
  int k = sign(b.z);
  if (cs[i,0,0] && cs[0,j,0] && cs[i,j,0] &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k] ) {
    double val_0, val_k;
    // bilinear interpolation in x-y-planes when all neighbors are defined
    val_0 = (s[0,0,0]*(1. - fabs(b.x)) + s[i,0,0]*fabs(b.x))*(1. - fabs(b.y)) +
      (s[0,j,0]*(1. - fabs(b.x)) + s[i,j,0]*fabs(b.x))*fabs(b.y);
    val_k = (s[0,0,k]*(1. - fabs(b.x)) + s[i,0,k]*fabs(b.x))*(1. - fabs(b.y)) +
      (s[0,j,k]*(1. - fabs(b.x)) + s[i,j,k]*fabs(b.x))*fabs(b.y);
    // trilinear interpolation when all neighbors are defined
    return (val_0*(1. - fabs(b.z)) + val_k*fabs(b.z));
  }
#endif 
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(b.x);
      if (cs[i] && cs[])
	val += fabs(b.x)*(s[i] - s[]);
      else if (cs[-i] && cs[] )
	val += fabs(b.x)*(s[] - s[-i]);
    }
    return val;
  }
}


#define CS_FLOOR 0.1

extern scalar * scalars;
extern vector * myvectors;

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

// #include "predictor-corrector.h"
#include "local-timestepping.h"

extern scalar rho, E;
extern vector w;
extern double gammao;

static void advance_lts_compressible (scalar * output, scalar * input,
                                      scalar * updates, double dt_global)
{
  foreach() {
    if (cs[] <= 0.) {
      scalar o;
      for (o in output)
        o[] = 0.;
    continue;
    } 
  
    double rho_v = rho[], w2 = 0.;
    foreach_dimension()
      w2 += sq(w.x[]);

    double dt_local = dt_global;
    if (rho_v > 1e-30) {
      double p_v = (gammao - 1.) * (E[] - 0.5 * w2 / rho_v);
      if (p_v < 0.) p_v = 1e-30;
      double a = sqrt(w2) / rho_v + sqrt(gammao * p_v / rho_v);
      if (a > 0.)
        dt_local = CFL * cm[] * Delta / a;
    }
    dt_local = fmin(dt_local, dt_global);

    scalar o, i, u;
    for (o,i,u in output, input, updates)
      o[] = i[] + dt_local * u[];
  }

  boundary (evolving);
}

/**
as well as a function which, given the state of each quantity,
returns the fluxes and the minimum/maximum eigenvalues
(i.e. `eigenvalue[0]`/`eigenvalue[1]`) of the corresponding Riemann
problem. */

void flux (const double * state, double * flux, double * eigenvalue);

/**
The generic time-integration scheme in `predictor-corrector.h` needs
to know which fields are updated i.e. all the scalars + the components
of all the vector fields. It also needs a function to compute the
time-derivatives of the evolving variables. */

scalar * evolving;
double update_conservation (scalar * conserved, scalar * updates, double dtmax);

event defaults (i = 0)
{
  evolving = list_concat (scalars, (scalar *) myvectors);
  update = update_conservation;
  lts_advance = advance_lts_compressible;

  /**
  We switch to a pure minmod limiter by default for increased
  robustness. */
  
  theta = 1.0;

  /**
  With the MUSCL scheme we use the CFL depends on the dimension of the
  problem. */

  if (CFL > 1./dimension)
    CFL = 1./dimension;
  
  /**
  On trees we need to replace the default bilinear
  refinement/prolongation with linear so that reconstructed values
  also use slope limiting. */

  #if TREE
  for (scalar s in evolving) {
    s.refine = s.prolongation = refine_embed_linear;
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
  #endif
}

/**
At the end of the run we need to free the list (to avoid a memory
leak). */

event cleanup (i = end) free (evolving);

/**
User initialisation happens here. */

event init (i = 0);

/**
### Computing fluxes

The core of the central-upwind scheme (see e.g. section 3.1 of
[Kurganov & Levy, 2002](references.bib#kurganov2002)) is the
approximate solution of the Riemann problem given by the left and
right states to get the fluxes `f`. */

void mirror_state(const double * fluid, double * ghost, const double * n, int len) {
  ghost[0] = fluid[0]; // rho
  ghost[1] = fluid[1]; // E

  double wdotn = 0.;
  for (int d = 0; d < dimension; d++)
      wdotn += fluid[2 + d] * n[d];
  for (int d = 0; d < dimension; d++)
      ghost[2 + d] = fluid[2 + d] - 2.*wdotn*n[d];
}

double riemann_Roe(const double * right, const double * left,
  double Delta, double * f, int len,
  double dtmax) {
  double rhoL = left[0], rhoR = right[0];
  double EL   = left[1], ER   = right[1];
  double wnL  = left[2], wnR  = right[2];

  if (rhoL < SEPS && rhoR < SEPS) {
      for (int i = 0; i < len; i++) f[i] = 0.;
      return dtmax;
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

    double a = fmax(fabs(lam[0]), fabs(lam[3]));
    if (a > 0.) {
        double dt = CFL * Delta / a;
        if (dt < dtmax)
            dtmax = dt;
    }
  return dtmax;         
}

static double riemann (const double * right, const double * left,
		       double Delta, double * f, int len, 
		       double dtmax)
{
  double fr[len], fl[len], er[2], el[2];
  flux (right, fr, er);
  flux (left,  fl, el);
  double ap = max(er[1], el[1]); ap = max(ap, 0.);
  double am = min(er[0], el[0]); am = min(am, 0.);
  double a = max(ap, -am); 

  if (a > 0.) {
    for (int i = 0; i < len; i++)
      f[i] = (ap*fl[i] - am*fr[i] + ap*am*(right[i] - left[i]))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < dtmax)
      dtmax = dt;
  }
  else
    for (int i = 0; i < len; i++)
      f[i] = 0.;
  return dtmax;
}


/**
## Utilities */

#define SEPS 1e-30                    

/**
## Computing fluxes and updates */

double update_conservation (scalar * conserved, scalar * updates, double dtmax)
{
  foreach() {
    if (cs[] <= 0.)
      for (scalar s in conserved)
        s[] = 0.;
  }

  vector * slopes = NULL;
  for (scalar s in conserved) {
    vector slope = new vector;
    foreach_dimension() {
      slope.x.gradient = zero;
    #if TREE
      slope.x.prolongation = refine_embed_linear;
    #endif
    }
    slopes = vectors_append (slopes, slope);
  }
  // gradients (conserved, slopes);
  foreach() {
    scalar s; vector g;
    for (s,g in conserved, slopes) {
      foreach_dimension() {
        if (!fs.x[] && !fs.x[1])
          g.x[] = 0.;
        else if (!fs.x[])
          g.x[] = s.gradient(s[], s[], s[1])/Delta;
        else if (!fs.x[1])
          g.x[] = s.gradient(s[-1], s[], s[])/Delta;
        else
          g.x[] = s.gradient(s[-1], s[], s[1])/Delta;
      }
    }
  }

  vector * lflux = NULL;
  int len = list_len (conserved);
  for (scalar s in conserved) {
    vector f1 = new face vector;
    lflux = vectors_append (lflux, f1);
  }

  int scalars_len = list_len (scalars);

  scalar * scalars = list_copy (conserved);
  if (scalars) scalars[scalars_len].i = -1;
  vector * vectors = vectors_from_scalars (&conserved[scalars_len]);
  
  vector * scalar_slopes = vectors_copy (slopes);
  if (scalar_slopes) scalar_slopes[scalars_len] = (vector){{-1}};
  tensor * vector_slopes = tensors_from_vectors (&slopes[scalars_len]);

  vector * scalar_fluxes = vectors_copy (lflux);
  if (scalar_fluxes) scalar_fluxes[scalars_len] = (vector){{-1}};
  tensor * vector_fluxes = tensors_from_vectors (&lflux[scalars_len]);

  for (vector f in lflux)
    foreach_face()
      f.x[] = 0.;

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

  foreach_face (reduction (min:dtmax)) {

    double r[len], l[len];
    double f[len];
    double dx = Delta/2.;
    int i = 0;
    scalar s;
    vector g;
    for (s,g in scalars,scalar_slopes) {
      r[i] = s[] - dx*g.x[];
      l[i++] = s[-1] + dx*g.x[-1];
    }
    vector v;
    tensor t;
    for (v,t in vectors,vector_slopes) {      
      r[i] = v.x[] - dx*t.x.x[];
      l[i++] = v.x[-1] + dx*t.x.x[-1];
    #if dimension > 1
        r[i] = v.y[] - dx*t.y.x[];
        l[i++] = v.y[-1] + dx*t.y.x[-1];
    #endif
    #if dimension > 2
        r[i] = v.z[] - dx*t.z.x[];
        l[i++] = v.z[-1] + dx*t.z.x[-1];
    #endif
    }

    double eff_delta = (fm.x[] > SEPS) ? Delta*cm[]/fm.x[] : Delta;

    bool right_valid = (cs[] > 0.);
    bool left_valid = (cs[-1] > 0.);

    if (!right_valid && !left_valid) {
      for (int k = 0; k < len; k++) f[k] = 0.;
    }
    else if (left_valid && !right_valid) {
      double n_wall[dimension];
      n_wall[0] = n_embed.x[];
#if dimension > 1
      n_wall[1] = n_embed.y[];
#endif
#if dimension > 2
      n_wall[2] = n_embed.z[];
#endif
      mirror_state(l, r, n_wall, len);
      riemann_Roe(r, l, eff_delta, f, len, dtmax);
    }
    else if (right_valid && !left_valid) {
      double n_wall[dimension];
      n_wall[0] = n_embed.x[-1];
#if dimension > 1
      n_wall[1] = n_embed.y[-1];
#endif
#if dimension > 2
      n_wall[2] = n_embed.z[-1];
#endif
      mirror_state(r, l, n_wall, len);
      riemann_Roe(r, l, eff_delta, f, len, dtmax);
    }
    else {
      double local_dtmax = riemann_Roe(r, l, eff_delta, f, len, dtmax);
      if (cs[] >= CS_FLOOR && cs[-1] >= CS_FLOOR)
        dtmax = fmin(dtmax, local_dtmax);
    }

    i = 0;
    for (vector fl in scalar_fluxes)
      fl.x[] = fm.x[]*f[i++];
    for (tensor fv in vector_fluxes) {
      fv.x.x[] = fm.x[]*f[i++];
      #if dimension > 1
        fv.y.x[] = fm.x[]*f[i++];
      #endif
      #if dimension > 2
        fv.z.x[] = fm.x[]*f[i++];
      #endif
    }
  }

  foreach() {
    if (cs[] <= 0.) {
      for (scalar ds in updates)
        ds[] = 0.;
      continue;
    }
    
    scalar ds;
    vector f;
    int idx = 0;
    for (ds,f in updates,lflux) {
      ds[] = 0.;
      foreach_dimension()
	      ds[] += (f.x[] - f.x[1])/(cm[]*Delta);
      idx++;
    }

    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry(point, &b, &n);
      
      double rho_b = my_embed_interpolate(point, rho, b);
      double E_b   = my_embed_interpolate(point, E, b);
      double wx_b  = my_embed_interpolate(point, w.x, b);
      double w2_b  = sq(wx_b);
    #if dimension > 1
      double wy_b  = my_embed_interpolate(point, w.y, b);
      w2_b += sq(wy_b);
    #endif
    #if dimension > 2
      double wz_b  = my_embed_interpolate(point, w.z, b);
      w2_b += sq(wz_b);
    #endif

    double p_val = 0.;
      if (rho_b > SEPS) {
        p_val = (gammao - 1.) * (E_b - 0.5 * w2_b / rho_b);
        if (p_val < 0. || p_val != p_val) p_val = 0.;
      } else {
        p_val = 0.;
      }

      // double c_f = sqrt(gammao * p_val / fmax(rho_b, SEPS));
      // double un_f = (wx_b * n.x + wy_b * n.y) / fmax(rho_b, SEPS);
      // // // 精确等熵 Riemann 不变量
      // double c_wall = c_f + 0.5 * (gammao - 1.) * un_f;
      // if (c_wall < SEPS) c_wall = c_f;
      // double p_wall = p_val * pow(c_wall / c_f, 2. * gammao / (gammao - 1.));
      // p_wall = fmax(p_wall, 0.);
      double p_wall = p_val;
      
      idx = 0;
      for (scalar ds in updates) {
        if (idx >= scalars_len && idx < scalars_len + dimension) {
          int dim_idx = idx - scalars_len;
          double n_comp = (dim_idx == 0) ? n.x :
          #if dimension > 1
                          (dim_idx == 1) ? n.y :
          #endif
          #if dimension > 2
                          (dim_idx == 2) ? n.z :
          #endif
                          0.;
          ds[] -= p_wall * area * n_comp / (cm[] * Delta);
        }
        idx++;
      }
    }
  }

  free (scalars);
  free (vectors);
  free (scalar_slopes);
  free (vector_slopes);
  free (scalar_fluxes);
  free (vector_fluxes);
  delete ((scalar *) slopes); free (slopes);
  delete ((scalar *) lflux); free (lflux);
  delete ((scalar *){n_embed});
  return dtmax;
}
