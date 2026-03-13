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

#define CS_FLOOR 0.1

extern scalar * scalars;
extern vector * myvectors;

extern scalar rho, E;
extern vector w;
extern double gammao;

/**
as well as a function which, given the state of each quantity,
returns the fluxes and the minimum/maximum eigenvalues
(i.e. `eigenvalue[0]`/`eigenvalue[1]`) of the corresponding Riemann
problem. */

void flux (const double * state, double * flux, double * eigenvalue);

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

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

  /**
  We switch to a pure minmod limiter by default for increased
  robustness. */
  
  theta = 1.;

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
  gradients (conserved, slopes);

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
    double local_dtmax = riemann (r, l, eff_delta, f, len, HUGE);
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
      
      double rho_b = embed_interpolate(point, rho, b);
      double E_b   = embed_interpolate(point, E, b);
      double wx_b  = embed_interpolate(point, w.x, b);
      double w2_b  = sq(wx_b);
    #if dimension > 1
      double wy_b  = embed_interpolate(point, w.y, b);
      w2_b += sq(wy_b);
    #endif
    #if dimension > 2
      double wz_b  = embed_interpolate(point, w.z, b);
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
          ds[] -= p_val * area * n_comp / (cm[] * Delta);
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
