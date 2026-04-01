/**
# A generic solver for systems of conservation laws (Ghost Cell Mirror version)

Based on Kurganov-Tadmor central-upwind scheme with ghost cell (mirror state)
treatment for embedded boundaries.
*/

double my_embed_interpolate (Point point, scalar s, coord b)
{
  int i = sign(b.x), j = sign(b.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j])
    return ((s[]*(1. - fabs(b.x)) + s[i]*fabs(b.x))*(1. - fabs(b.y)) + 
	    (s[0,j]*(1. - fabs(b.x)) + s[i,j]*fabs(b.x))*fabs(b.y));
#else // dimension == 3
  int k = sign(b.z);
  if (cs[i,0,0] && cs[0,j,0] && cs[i,j,0] &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k] ) {
    double val_0, val_k;
    val_0 = (s[0,0,0]*(1. - fabs(b.x)) + s[i,0,0]*fabs(b.x))*(1. - fabs(b.y)) +
      (s[0,j,0]*(1. - fabs(b.x)) + s[i,j,0]*fabs(b.x))*fabs(b.y);
    val_k = (s[0,0,k]*(1. - fabs(b.x)) + s[i,0,k]*fabs(b.x))*(1. - fabs(b.y)) +
      (s[0,j,k]*(1. - fabs(b.x)) + s[i,j,k]*fabs(b.x))*fabs(b.y);
    return (val_0*(1. - fabs(b.z)) + val_k*fabs(b.z));
  }
#endif 
  else {
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

extern scalar rho, E;
extern vector w;
extern double gammao;

void flux (const double * state, double * flux, double * eigenvalue);

#include "predictor-corrector.h"

scalar * evolving;
double update_conservation (scalar * conserved, scalar * updates, double dtmax);

event defaults (i = 0)
{
  evolving = list_concat (scalars, (scalar *) myvectors);
  update = update_conservation;
  theta = 1.;

  if (CFL > 1./dimension)
    CFL = 1./dimension;
  
#if TREE
  for (scalar s in evolving) {
    s.refine = s.prolongation = refine_embed_linear;
    s.restriction = restriction_volume_average;
    s.dirty = true;
  }
#endif
}

event cleanup (i = end) free (evolving);

event init (i = 0);

/**
## Riemann solver */

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

#define SEPS 1e-30                    

/**
## Construct mirror (ghost) state for slip wall

The normal is passed as a plain double array (NOT coord) so that
foreach_face dimension rotation keeps components consistent with
the state array ordering.

  w_ghost = w_fluid - 2*(w_fluid . n)*n
*/

static void mirror_state (const double * fluid, double * ghost,
                          const double * n, int len, int scalars_len)
{
  for (int i = 0; i < scalars_len; i++)
    ghost[i] = fluid[i];

  double wdotn = 0.;
  for (int d = 0; d < dimension; d++)
    wdotn += fluid[scalars_len + d] * n[d];

  for (int d = 0; d < dimension; d++)
    ghost[scalars_len + d] = fluid[scalars_len + d] - 2. * wdotn * n[d];
}

/**
## Computing fluxes and updates */

double update_conservation (scalar * conserved, scalar * updates, double dtmax)
{
  foreach() {
    if (cs[] <= 0.)
      for (scalar s in conserved)
        s[] = 0.;
  }

  /**
  Precompute embed normals for all cut cells into a vector field,
  so that they can be accessed with index offsets in foreach_face. */

  vector n_embed = new vector;
  foreach() {
    foreach_dimension()
      n_embed.x[] = 0.;
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      embed_geometry (point, &b, &n);
      foreach_dimension()
        n_embed.x[] = n.x;
    }
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

    bool right_valid = (cs[] >= CS_FLOOR);
    bool left_valid  = (cs[-1] >= CS_FLOOR);

    if (!right_valid && !left_valid) {
      for (int k = 0; k < len; k++)
        f[k] = 0.;
      i = 0;
      for (vector fl in scalar_fluxes)
        fl.x[] = fm.x[] * f[i++];
      for (tensor fv in vector_fluxes) {
        fv.x.x[] = fm.x[] * f[i++];
      #if dimension > 1
        fv.y.x[] = fm.x[] * f[i++];
      #endif
      #if dimension > 2
        fv.z.x[] = fm.x[] * f[i++];
      #endif
      }
      continue;
    }

    /**
    Standard left/right MUSCL reconstruction. */
    scalar s;
    vector g;
    i = 0;
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

    if (!right_valid && left_valid) {
      /**
      Right side is solid => mirror the left (fluid) state.
      Use a plain double array for the normal so that dimension
      rotation keeps it consistent with the state vector ordering. */
      double n_wall[dimension];
      double nmag2 = 0.;

      n_wall[0] = n_embed.x[];
#if dimension > 1
      n_wall[1] = n_embed.y[];
#endif
#if dimension > 2
      n_wall[2] = n_embed.z[];
#endif
      for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);

      if (nmag2 < SEPS) {
        n_wall[0] = n_embed.x[-1];
#if dimension > 1
        n_wall[1] = n_embed.y[-1];
#endif
#if dimension > 2
        n_wall[2] = n_embed.z[-1];
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
      mirror_state (l, r, n_wall, len, scalars_len);
    }
    else if (right_valid && !left_valid) {
      /**
      Left side is solid => mirror the right (fluid) state. */
      double n_wall[dimension];
      double nmag2 = 0.;

      n_wall[0] = n_embed.x[-1];
#if dimension > 1
      n_wall[1] = n_embed.y[-1];
#endif
#if dimension > 2
      n_wall[2] = n_embed.z[-1];
#endif
      for (int d = 0; d < dimension; d++) nmag2 += sq(n_wall[d]);

      if (nmag2 < SEPS) {
        n_wall[0] = n_embed.x[];
#if dimension > 1
        n_wall[1] = n_embed.y[];
#endif
#if dimension > 2
        n_wall[2] = n_embed.z[];
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
      mirror_state (r, l, n_wall, len, scalars_len);
    }

    double eff_delta = (fm.x[] > SEPS) ? Delta*cm[]/fm.x[] : Delta;
    double local_dtmax = riemann (r, l, eff_delta, f, len, HUGE);
    if (right_valid && left_valid)
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
    for (ds,f in updates,lflux) {
      ds[] = 0.;
      foreach_dimension()
	ds[] += (f.x[] - f.x[1])/(cm[]*Delta);
    }
  }

  delete ((scalar *){n_embed});
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

event cell_linking (i++) {
  for (scalar s in evolving) {
    foreach() {
      if (cs[] > 0. && cs[] < CS_FLOOR) {
        double sum_val = 0., sum_w = 0.;
        foreach_neighbor(1) {
          if (cs[] >= CS_FLOOR) {
            sum_val += cs[] * s[];
            sum_w += cs[];
          }
        }
        if (sum_w > 0.)
          s[] = sum_val / sum_w;
      }
    }
  }
}