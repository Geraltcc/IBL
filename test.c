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
#define ibl_alpha 0.2
#define ibl_CFL 0.1
#define ibl_niter 100000
// #define ibl_niter 0

#define TE_RADIUS 2.0

scalar d[];
coord* p = NULL;

scalar * dlist = {rho, w.x, w.y, E, cs};

scalar ue[];
scalar ibl_theta[];
vector ue_vec[];
vector grad_ue[];
vector grad_theta[];
scalar Cf2[];
scalar Me[];

void solve_ibl();
double abs_max(scalar s);
double s_max(scalar s);
double s_min(scalar s);
void plot_vof(char* name, char *filename);
double turbulent_Cf2(double Hk, double Re_th, double Me);


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

    build_cutcell_cache();

    int count = 0;
    foreach_cache(cutcells, reduction(+:count)) {
        int cnt = 0;
        foreach_neighbor(1) {
            if (is_solid) cnt++;
        }
        if (cnt == 0) count++;
    }
    fprintf(stderr, "nested cutcells count: %d\n", count);

    // ========== 1. 计算 u_e 和流向 ==========
    foreach_cache(cutcells) {
        coord n, b;
        embed_geometry(point, &b, &n);
        double rho_s = my_embed_interpolate(point, rho, b);
        double ux = my_embed_interpolate(point, w.x, b) / rho_s;
        double uy = my_embed_interpolate(point, w.y, b) / rho_s;
        double un = ux * n.x + uy * n.y;
        double etx = ux - un * n.x;
        double ety = uy - un * n.y;
        ue[] = sqrt(etx * etx + ety * ety);
        ue_vec.x[] = etx;
        ue_vec.y[] = ety;

        double E_s = my_embed_interpolate(point, E, b);
        double p_s = (gammao - 1.) * (E_s - 0.5*(sq(ux) + sq(uy))/rho_s);
        double a_s = sqrt(gammao * p_s / rho_s);
        Me[] = ue[] / fmax(a_s, 1e-10);
    }

    // foreach_cache(cutcells) {
    //     double ux = w.x[] / rho[];
    //     double uy = w.y[] / rho[];
    //     ue[] = sqrt(sq(ux) + sq(uy));
    // }
    boundary({ue, ue_vec});

    FILE* f_ue = fopen("data/ue.dat", "w");
    fprintf(f_ue, "# x y ue\n");
    foreach_cache(cutcells) {
        // double dist = sqrt(sq(x - stag_x) + sq(y - stag_y));

        fprintf(f_ue, "%.6e %.6e\n", x, w.x[]);
    }
    fflush(f_ue);
    fclose(f_ue);
    system("gnuplot plot_ue.gp");

    double Delta_min = L0 / (1 << maxlevel);
    double dtau = ibl_CFL * Delta_min;
    double te_radius = TE_RADIUS * Delta_min;

    // ========== 3. 驻点识别 ==========
    double stag_x = 0., stag_y = 0., ue_min = HUGE;
    foreach_cache(cutcells) {
        if (ue[] < ue_min) {
            ue_min = ue[];
            stag_x = x;
            stag_y = y;
        }
    }

    double K_stag = 0.;
    int k_count = 0.;
    foreach_cache(cutcells) {
        double dist = sqrt(sq(x - stag_x) + sq(y - stag_y));
        if (dist > 0.5*Delta && dist < 3.*Delta && ue[] > 1e-6) {
            K_stag += ue[] / dist;
            k_count++;
        }
    }
    if (k_count > 0) K_stag /= k_count;
    else K_stag = 1.;

    double Cf2_stag = turbulent_Cf2(1.3, 100., 0.);
    double theta_stag = sqrt(2.*Cf2_stag*ibl_nu / ((2. + ibl_H)*K_stag));
    fprintf(stderr, "IBL: stag=(%.4f, %.4f) ue_min=%.4e K=%.2f theta_stag=%.4e\n", 
        stag_x, stag_y, ue_min, K_stag, theta_stag);

    // 尾缘识别
    double te_x = stag_x, te_y = stag_y, dist_max = 0.;
    foreach_cache(cutcells) {
        double d = sq(x - stag_x) + sq(y - stag_y);
        if (d >dist_max) {
            dist_max = d;
            te_x = x;
            te_y = y;
        }
    }

    fprintf(stderr, "IBL: te=(%.4f, %.4f)\n", te_x, te_y);

    
#define smooth_grad(g) do { \
    foreach_cache(cutcells) { \
        double dist_te = sqrt(sq(x - te_x) + sq(y - te_y)); \
        if (dist_te < te_radius) { \
            double w = 0.25 * (1.0 - cos(M_PI *dist_te / te_radius)); \
            foreach_dimension() \
                g.x[] *= w; \
        }       \
    } \
} while (0)

#define clean_grad(g) do { \
    foreach_cache(cutcells) { \
        double dist_te = sqrt(sq(x - te_x) + sq(y - te_y)); \
        if (dist_te < 3.*Delta_min) { \
            foreach_dimension() \
                g.x[] = 0.; \
        }       \
    } \
} while (0)

#define clean_scalar(s) do { \
    foreach_cache(cutcells) { \
        double dist_te = sqrt(sq(x - te_x) + sq(y - te_y)); \
        if (dist_te < 3.0*Delta_min) { \
            s[] = 0.; \
        } \
    } \
} while (0)

    // ========== 2. 计算 ∇_s u_e (Euler 解固定，只算一次) ==========
    cutcell_tangential_gradient_lsq_backward(ue, grad_ue);
    // cutcell_tangential_gradient_lsq_forward_fixed(ue, grad_ue);
    // smooth_grad(grad_ue);
    // clean_grad(grad_ue);

    // ========== 4. 初始化 θ ==========

    foreach_cache(cutcells)
        ibl_theta[] = theta_stag;
    boundary({ibl_theta}); 

    // ========== 5. 伪时间迭代 ==========

    fprintf(stderr, "IBL: Delta_min=%.4e dtau=%.4e niter=%d\n",
        Delta_min, dtau, ibl_niter);

    double ibl_theta_min, ibl_theta_max;

    for (int iter = 0; iter < ibl_niter; iter++) {
        cutcell_tangential_gradient_lsq_backward(ibl_theta, grad_theta);
        // cutcell_tangential_gradient_lsq_forward_fixed(ibl_theta, grad_theta);
        // smooth_grad(grad_theta);
        // clean_grad(grad_theta);

        double res_max = 0.;
        foreach_cache(cutcells) {
            double dist = sqrt(sq(x - stag_x) + sq(y - stag_y));
            if (dist < 1.5*Delta) {
                ibl_theta[] = theta_stag;
                continue;
            }
    
            double emag = sqrt(sq(ue_vec.x[]) + sq(ue_vec.y[]));
            if (emag < 1e-10) continue;
            double ex = ue_vec.x[] / emag;
            double ey = ue_vec.y[] / emag;

            // double theta_here = ibl_theta[];
            // double theta_up = theta_here;
            // double dist_up = 0.;
            // double dot_min = -1e-10;

            // double x0 = x, y0 = y;

            // foreach_neighbor(1) {
            //     if (is_cutcell) {
            //         double dx = x - x0, dy = y - y0;
            //         double dot = dx*ex + dy*ey;
            //         double dist = sqrt(sq(dx) + sq(dy));
            //         if (dot < dot_min && dist > 1e-10) {
            //             dot_min = dot;
            //             theta_up = ibl_theta[];
            //             dist_up = dist;
            //         }
            //     }
            // }

            // double dtheta_dxi = (dist_up > 1e-10) ? 
            //     (theta_here - theta_up) / dist_up : 0.;


            double dtheta_dxi = ex * grad_theta.x[] + ey * grad_theta.y[];
            double due_dxi    = ex * grad_ue.x[]    + ey * grad_ue.y[];

            double ue_safe = fmax(ue[], 1e-8);
            double Re_th = ue_safe * ibl_theta[] / ibl_nu;
            double Hk = (ibl_H - 0.290 * sq(Me[])) / (1. + 0.113 * sq(Me[]));
            Cf2[] = turbulent_Cf2(Hk, Re_th, Me[]);

            double sum_th = 0.;
            int cnt = 0;
            
            double x0 = x, y0 = y;
            foreach_neighbor(1) {
                double dx = x - x0, dy = y - y0;
                double dot = dx*ex + dy*ey;
                if (is_cutcell && dot < 0.) {
                    sum_th += ibl_theta[];
                    cnt++;
                }
            }

            double th_avg = (cnt > 0) ? 
                (1. - ibl_alpha)*ibl_theta[] + ibl_alpha*sum_th/cnt : ibl_theta[];

            double theta_new = th_avg + dtau * (Cf2[] - dtheta_dxi - (2. + ibl_H) * ibl_theta[] / ue_safe * due_dxi);

            // double theta_new = ibl_theta[] + dtau * (Cf2[] - dtheta_dxi - (2. + ibl_H) * ibl_theta[] / ue_safe * due_dxi);
            if (theta_new < 1e-12) theta_new = 1e-12;

            double res = fabs(theta_new - ibl_theta[]);
            if (res > res_max) res_max = res;

            // double dist_te = sqrt(sq(x - te_x) + sq(y - te_y));
            // if (dist_te > te_radius)
            ibl_theta[] = theta_new;
        }

        // clean_scalar(ibl_theta);

        if (iter % 1000 == 0) {
            fprintf(stderr, "[%05d] res_max=%.4e\n", iter, res_max);
            // if (res_max < 1e-10)
            //     break;

            FILE* f_theta = fopen("data/theta.dat", "w");
            fprintf(f_theta, "# x y theta\n");
            foreach_cache(cutcells) {
                if (y > 0.) {
                    double dist = sqrt(sq(x - stag_x) + sq(y - stag_y));

                    fprintf(f_theta, "%.6e %.6e\n", dist, ibl_theta[]);
                }
            }
            fflush(f_theta);
            fclose(f_theta);

            FILE* f_data = fopen("data/Cf2.dat", "w");
            fprintf(f_data, "# x y Cf\n");
            foreach_cache(cutcells) {
                // double dist = sqrt(sq(x - stag_x) + sq(y - stag_y));
                if (y > 0.) {
                    double dist = fabs(x - stag_x);
            
                    fprintf(f_data, "%.6e %.6e\n", dist, Cf2[]*2.);
                }
            }
            fflush(f_data);
            fclose(f_data);

            system("gnuplot plot_th.gp");
            system("gnuplot plot_Cf.gp");
        }
    }
    boundary({ibl_theta});

    ibl_theta_min = s_min(ibl_theta);
    ibl_theta_max = s_max(ibl_theta);
    // view(width = 2000, height = 2000);
    // clear();
    //  ("cs", "fs", color = "ibl_theta", map=viridis, min = ibl_theta_min, max = ibl_theta_max, linear=false, lc = None);
    // colorbar(map = viridis, min = ibl_theta_min, max = ibl_theta_max);
    // save("output/ibl_theta.png");
    fprintf(stderr, "[info] theta_min=%.4e theta_max=%.4e\n", ibl_theta_min, ibl_theta_max);



    plot_vof("grad_theta.x", "output/grad_theta.png");
}

double turbulent_Cf2(double Hk, double Re_th, double Me) {
    double Fc = sqrt(1. + 0.2 * sq(Me));
    double Re_Fc = Re_th / Fc;
    if (Re_Fc < 10.) return 0.;
    double log_Re = log10(Re_Fc);
    double Cf_Fc = 0.3 * exp(-1.33 * Hk)
                 * pow(log_Re, -1.74 - 0.31 * Hk)
                 + 0.00011 * (tanh(4. - Hk / 0.875) - 1.);
    return fmax(0.5 * Cf_Fc / Fc, 0.);
}

double abs_max(scalar s) {
    double max = 0.;
    foreach(reduction(max:max)) {
        if (max < fabs(s[]))
            max = fabs(s[]);
    }
    return max;
}

double s_max(scalar s) {
    double max = -HUGE;
    foreach_cache(cutcells, reduction(max:max)) {
        if (max < s[])
            max = s[];
    }
    return max;
}

double s_min(scalar s) {
    double min = HUGE;
    foreach_cache(cutcells, reduction(min:min)) {
        if (min > s[])
            min = s[];
    }
    return min;
}

void plot_vof(char* name, char *filename) {
    scalar s = lookup_field (name);
    double m = abs_max(s);
    view(width = 2000, height = 2000);
    clear();
    squares(name, map = blue_white_red, min = -m, max = m, linear = false);
    colorbar(map = blue_white_red, min = -m, max = m);
    save(filename);
}
