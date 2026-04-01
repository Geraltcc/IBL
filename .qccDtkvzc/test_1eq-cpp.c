@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "test_1eq.c"
#include "grid/quadtree.h"
#include "embed.h"
#include "distance.h"
#include "mycompressible_nofr.h"
#include "view.h"
#include "cmap.h"
#include "ibl.h"

#define maxlevel 10
#define minlevel 6

#define mach 0.6
#define rho0 1.
#define p0 1.

#define ibl_nu 1e-7
#define ibl_H  1.5
#define ibl_ReCf2 0.221
#define ibl_alpha 0.5
#define ibl_CFL 0.1
#define ibl_niter 100000

scalar d[];
coord* p = NULL;

scalar * dlist = {rho, w.x, w.y, E, cs};

scalar ue[];
scalar ibl_theta[];
vector ue_vec[]
vector grad_ue[];
vector grad_theta[];
scalar Cf2[];
scalar Me[];

int main() {
    L0 = 5.0;
    size(L0);
    origin(-L0/2., -L0/2.);

    init_grid(1 << minlevel);

    restore("restart");
    boundary(all);
    
    FILE* fp = fopen("NACA0012_blunt.gnu", "r");
    p = input_xy(fp);
    distance(d, p);
    
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);
    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    foreach() {
        ue[] = 0.;
        ibl_theta[] = 0.;
        Cf2[] = 0.;
        foreach_dimension() {
            ue_vec.x[] = 0.;
            grad_ue.x[] = 0.;
            grad_theta.x[] = 0.;
        }
    }
}
#endif
