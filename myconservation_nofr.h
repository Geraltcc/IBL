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

#include "predictor-corrector.h"
// #include "local-timestepping.h"

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

    if (cs[] < CS_FLOOR)
      continue;
  
    double rho_v = rho[], w2 = 0.;
    foreach_dimension()
      w2 += sq(w.x[]);

    double dt_local = dt_global;
    if (rho_v > 1e-30) {
      double p_v = (gammao - 1.) * (E[] - 0.5 * w2 / rho_v);
      if (p_v < 0.) p_v = 1e-30;
      double a = sqrt(w2) / rho_v + sqrt(gammao * p_v / rho_v);
      if (a > 0.)
        dt_local = CFL * Delta / a;
    }

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
  // lts_advance = advance_lts_compressible;

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

static double riemann_hllc (double * right, double * left,
                             double Delta, double * f, int len,
                             double dtmax)
{
  // 两侧都无效 → 零通量
  double rhoL = left[0], rhoR = right[0];
  if ((rhoL <= 0. || rhoL != rhoL) && (rhoR <= 0. || rhoR != rhoR)) {
    for (int i = 0; i < len; i++)
      f[i] = 0.;
    return dtmax;
  }

  // 只有左侧无效 → 用右侧状态做零阶外推
  if (rhoL <= 0. || rhoL != rhoL) {
    for (int i = 0; i < len; i++)
      left[i] = right[i];   // 注意：left 需要是非 const 的，或者用局部拷贝
  }

  // 只有右侧无效 → 用左侧状态做零阶外推
  if (rhoR <= 0. || rhoR != rhoR) {
    for (int i = 0; i < len; i++)
      right[i] = left[i];
  }

  rhoL = left[0];
  rhoR = right[0];

  // 左状态
  double EL = left[1], wnL = left[2];
  double w2L = 0.;
  for (int i = 2; i < 2 + dimension; i++) w2L += sq(left[i]);
  double uL = wnL/rhoL;
  double pL = (gammao - 1.)*(EL - 0.5*w2L/rhoL);
  if (pL < 0.) pL = 0.;
  double cL = sqrt(gammao*pL/rhoL);

  // 右状态
  double ER = right[1], wnR = right[2];
  double w2R = 0.;
  for (int i = 2; i < 2 + dimension; i++) w2R += sq(right[i]);
  double uR = wnR/rhoR;
  double pR = (gammao - 1.)*(ER - 0.5*w2R/rhoR);
  if (pR < 0.) pR = 0.;
  double cR = sqrt(gammao*pR/rhoR);

  // 波速估计 (Einfeldt/Davis)
  double SL = min(uL - cL, uR - cR);
  double SR = max(uL + cL, uR + cR);

  // 接触波速 S*
  double Sstar = (pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR))
               / (rhoL*(SL - uL) - rhoR*(SR - uR));

  if (SL >= 0.) {
    // 全部取左通量
    f[0] = wnL;
    f[1] = uL*(EL + pL);
    f[2] = uL*wnL + pL;
    for (int i = 3; i < 2 + dimension; i++)
      f[i] = uL*left[i];
  }
  else if (SR <= 0.) {
    // 全部取右通量
    f[0] = wnR;
    f[1] = uR*(ER + pR);
    f[2] = uR*wnR + pR;
    for (int i = 3; i < 2 + dimension; i++)
      f[i] = uR*right[i];
  }
  else {
    // 中间状态
    double fL[len], fR[len];
    fL[0] = wnL;
    fL[1] = uL*(EL + pL);
    fL[2] = uL*wnL + pL;
    for (int i = 3; i < 2 + dimension; i++)
      fL[i] = uL*left[i];
    fR[0] = wnR;
    fR[1] = uR*(ER + pR);
    fR[2] = uR*wnR + pR;
    for (int i = 3; i < 2 + dimension; i++)
      fR[i] = uR*right[i];

    if (Sstar >= 0.) {
      // U*L 区域
      double coef = rhoL*(SL - uL)/(SL - Sstar);
      double UstarL[len];
      UstarL[0] = coef;
      UstarL[1] = coef*(EL/rhoL + (Sstar - uL)*(Sstar + pL/(rhoL*(SL - uL))));
      UstarL[2] = coef*Sstar;
      for (int i = 3; i < 2 + dimension; i++)
        UstarL[i] = coef*left[i]/rhoL;  // 切向动量守恒

      for (int i = 0; i < len; i++)
        f[i] = fL[i] + SL*(UstarL[i] - left[i]);
    }
    else {
      // U*R 区域
      double coef = rhoR*(SR - uR)/(SR - Sstar);
      double UstarR[len];
      UstarR[0] = coef;
      UstarR[1] = coef*(ER/rhoR + (Sstar - uR)*(Sstar + pR/(rhoR*(SR - uR))));
      UstarR[2] = coef*Sstar;
      for (int i = 3; i < 2 + dimension; i++)
        UstarR[i] = coef*right[i]/rhoR;

      for (int i = 0; i < len; i++)
        f[i] = fR[i] + SR*(UstarR[i] - right[i]);
    }
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > 0.) {
    double dt = CFL*Delta/a;
    if (dt < dtmax) dtmax = dt;
  }

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
    double local_dtmax = riemann_hllc (r, l, eff_delta, f, len, HUGE);
    // double local_dtmax = riemann (r, l, eff_delta, f, len, HUGE);
    if (cs[] >= CS_FLOOR && cs[-1] >= CS_FLOOR)
      dtmax = min(dtmax, local_dtmax);

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
    if (cs[] < CS_FLOOR) {
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
      }
      
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
          double cm_safe = fmax(cs[], CS_FLOOR);
          ds[] -= p_val * area * n_comp / (cm_safe * Delta);
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
  delete ((scalar *) slopes);
  free (slopes);
  delete ((scalar *) lflux);
  free (lflux);
  
  return dtmax;
}
