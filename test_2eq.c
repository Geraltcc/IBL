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
// #define ibl_H  1.5
#define ibl_ReCf2 0.221
#define ibl_alpha 0.5
#define ibl_CFL 0.1
#define ibl_dt 1e-4
#define ibl_niter 100

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
scalar ibl_H[];
// scalar H_star[];

scalar delta[];
vector grad_delta[];


int main() {
    L0 = 5.0;
    size(L0);
    origin(-L0/2., -L0/2.);

    init_grid(1 << minlevel);

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

    run();
}

rho[left] = dirichlet(rho0);
E[left] = dirichlet(p0/(gammao - 1.) + 0.5*rho0*sq(mach));
w.x[left] = dirichlet(mach*rho0);
w.y[left] = dirichlet(0.);

// 流出边界
rho[right] = neumann(0.);
E[right] = neumann(0.);
w.x[right] = neumann(0.);
w.y[right] = neumann(0.);

event init(i = 0) {
    // restore("restart");
    // boundary(all);
    
    FILE* fp = fopen("NACA0012_blunt.gnu", "r");
    p = input_xy(fp);
    distance(d, p);
    
    while(adapt_wavelet({d}, (double[]){1e-30}, maxlevel).nf);
    do {
        solid(cs, fs, (d[]+d[-1]+d[0,-1]+d[-1,-1])/4.);
    } while (adapt_wavelet({cs}, (double[]){1e-30}, maxlevel).nf);

    build_cutcell_cache();

    foreach() {
        if (cs[] > 0.) {
            rho[] = rho0;
            w.x[] = mach * rho0;
            w.y[] = 0.;
            E[] = p0/(gammao - 1.) + 0.5*rho0*sq(mach);
        } else {
            rho[] = 0.;
            w.x[] = 0.;
            w.y[] = 0.;
            E[] = 0.;
        }
    }
}

event adapt(i++) {
    adapt_wavelet({cs, rho, w}, (double[]){1e-30, 1e-2, 1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
    fractions_cleanup(cs, fs);
    build_cutcell_cache();
}

double w_max = 0.;
double w_min = 0.;
double rho_max = 0.;
double rho_min = 0.;
double cs_min = 0.;

double scalar_max(scalar s) {
    double max = -HUGE;
    foreach(reduction(max:max)) {
        if (max < s[]) {
            max = s[];
        }
    }
    return max;
}

double scalar_min(scalar s) {
    double min = HUGE;
    foreach(reduction(min:min)) {
        if (min > s[]) {
            min = s[];
        }
    }
    return min;
}

scalar wn[];

event logfile(i += 100) {
    double dw = 0.;
    if (dt > 1e-10)
        dw = change (w.x, wn)/dt;
    else
        dw = 0.;
    w_max = scalar_max(w.x);
    cs_min = HUGE;
    foreach(reduction(min:cs_min)) {
        if (cs[] > 0. && cs[] < cs_min) {
            cs_min = cs[];
        }
    }
    fprintf(stderr, "[log] i: %04d t: %.4e dt: %.4e grid: %ld w_max: %.4e cs_min: %.4e dw: %.4e\n", i, t, dt, grid->n, w_max, cs_min, dw);
}

double cal_Hk(double Ma, double H) {
    double Hk = (H - 0.290 * sq(Ma)) / (1. + 0.113 * sq(Ma));
    return fmax(Hk, 0.);
}

double cal_Cf2(double Ma, double Re, double Hk) {
    double Fc = sqrt(1. + 0.2*sq(Ma));
    return fmax(0.5 * 0.3 * exp(-1.33 * Hk) 
        * pow(log10(fmax(Re/Fc, 2.0)), -1.74 - 0.31*Hk)
        + 0.00011 * (tanh(4. - Hk/0.875) - 1.), 0.);
}

double cal_H_star(double Re, double Hk) {
    double H0 = (Re > 400.) ? (3. + 400. / Re) : 4.;
    double H_star = 0.;
    if (Hk < H0) {
        H_star = 1.505 + 4./Re + (0.165 - 1.6/sqrt(Re)) * pow(H0 - Hk, 1.6) / Hk;
    } else {
        H_star = 1.505 + 4./Re + sq(Hk - H0) * (0.04/Hk + 0.007 * log(Re)/sq(Hk - H0 + 4./log(Re)));
    }
    return fmax(H_star, 0.);
}

double cal_H_2star(double Ma, double Hk) {
    double H_2star = (0.064 / (Hk - 0.8) + 0.251) * sq(Ma);
    return fmax(H_2star, 0.);
}

double cal_Us(double H, double Hk, double H_star) {
    double Us = 0.5 * H_star * (1. - 4./3. * (Hk - 1)/H);
    return fmax(Us, 0.);
}

double cal_Cd(double Cf2, double Us, double Ctau) {
    double Cd = Cf2 * Us + Ctau * (1. - Us);
    return fmax(Cd, 0.);
}

double cal_Ctau_eq(double H, double Hk, double H_star, double Us) {
    double Ctau_eq = H_star * 0.015 / (1. - Us) * sq(Hk - 1.) / sq(Hk) / H;
    return fmax(Ctau_eq, 0.); 
}

void update_ue() {
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
        // ue_vec.x[] = etx;
        // ue_vec.y[] = ety;
        double norm = sqrt(sq(etx) + sq(ety));
        ue_vec.x[] = etx / norm;
        ue_vec.y[] = ety / norm;

        double E_s = my_embed_interpolate(point, E, b);
        double p_s = (gammao - 1.) * (E_s - 0.5*(sq(ux) + sq(uy))/rho_s);
        double a_s = sqrt(gammao * p_s / rho_s);
        Me[] = ue[] / fmax(a_s, 1e-10);
    }
    boundary({ue, ue_vec, Me});
}

double stag_x = 0.;
double stag_y = 0.;
double ue_min = HUGE;
double tail_x = 0.;

event ibl_theta_init(i = 0) {
    // get ue
    update_ue();

    foreach_cache(cutcells) {
        if (ue[] < ue_min) {
            ue_min = ue[];
            stag_x = x;
            stag_y = y;
        }
        ibl_theta[] = 1e-6;
        ibl_H[] = 1.4;
    }

    tail_x = -HUGE;
    foreach_cache(cutcells) { 
        if (x > tail_x) {
            tail_x = x;
        }   
    }

    for (scalar s in {ue, ibl_theta, Me, delta, ue_vec.x, ue_vec.y}) {
        s.refine = s.prolongation = refine_embed_linear;
        s.restriction = restriction_volume_average;
    }
}

double upwind_grad(scalar s, coord e, Point point) {
    coord x0 = {x, y};

    double s0 = s[];
    double sum = 0., w = 0.;

    foreach_neighbor(1) {
        if (!is_cutcell) continue;
        
        coord dx = {x - x0.x, y - x0.y};
        double dist2 = sq(dx.x) + sq(dx.y);
        
        if (dist2 < 1e-10) continue;

        double dot = dx.x * e.x + dx.y * e.y;

        if (dot < 0.) {
            double ds = -dot;
            double wt = ds / dist2;
            double slope = (s0 - s[]) / ds;
            sum += wt * slope;
            w += wt;
        }
    }
    return (w > 0.) ? sum / w : 0.;
}

double central_grad (scalar s, coord e, Point point) {
    coord x0 = {x, y};
    
    double s0 = s[];
    double sum = 0., wt_sum = 0.;

    foreach_neighbor(1) {
        if (!is_cutcell) continue;

        coord dx = {x - x0.x, y - x0.y};

        double dist2 = sq(dx.x) + sq(dx.y);
        if (dist2 < 1e-10) continue;

        double dot = dx.x * e.x + dx.y * e.y;

        if (fabs(dot) > 1e-20) {
            double wt = 1. / dist2;
            sum += wt * (s[] - s0) / dot;
            wt_sum += wt;
        }
    }
    return (wt_sum > 0.) ? sum / wt_sum : 0.;
}

#define clamp(val, edge) val = fmax(fmin(val, edge), -edge)

#define front_tail_handling() \
    double dist_stag = sqrt(sq(x - stag_x) + sq(y - stag_y)); \
    double dist_tail = fabs(x - tail_x); \
    if (dist_stag < 3.0 * Delta || dist_tail < 3.0 * Delta) { \
        ibl_theta[] = 1e-6; \
        continue; \
    } \


event ibl_theta_equation(i++) {
    update_ue();

    if (t > 0.2) {

    // 更新 Hstar
    for (int iter = 0; iter < ibl_niter; iter++) {

        // 计算 theta 方程
        scalar theta_tmp[];
        scalar H_tmp[];
        foreach_cache(cutcells) {
            front_tail_handling();

            // double ue_safe = fmax(ue[], 1e-8);
            // double Re_th = ue_safe * ibl_theta[] / ibl_nu;
            // double Hk = cal_Hk(Me[], H[]);

            coord e = {ue_vec.x[], ue_vec.y[]};
            double dtheta_dxi = upwind_grad(ibl_theta, e, point);
            double due_dxi = central_grad(ue, e, point);
            double ue_safe = fmax(ue[], 1e-8);
            
            double theta_new = ibl_theta[] + ibl_dt * (
                Cf2[] - dtheta_dxi - 
                (2. + ibl_H[] - sq(Me[])) * ibl_theta[] / ue_safe * due_dxi
            );
            theta_new = fmax(theta_new, 0.);
            theta_tmp[] = theta_new;

            // ================ H Equation =================
            double Re_th = fmax(rho[] * ue_safe * ibl_theta[] / ibl_nu, 1e-4);
            double Hk = fmax(cal_Hk(Me[], ibl_H[]), 1e-4);
            double H_star = cal_H_star(Re_th, Hk);
            Cf2[] = cal_Cf2(Me[], Re_th, Hk);
            double Us = cal_Us(ibl_H[], Hk, H_star);
            double Ctau = cal_Ctau_eq(ibl_H[], Hk, H_star, Us);
            double Cd = cal_Cd(Cf2[], Us, Ctau);

            double H_2star = cal_H_2star(Me[], Hk);

            // double dH_dxi = upwind_grad(ibl_H, e, point);

            double eps = 1e-6;
            double H_p = ibl_H[] + eps;
            double Hk_p = fmax(cal_Hk(Me[], H_p), 1e-4);
            double H_star_p = cal_H_star(Re_th, Hk_p);
            double dHstar_dH = (H_star_p - H_star) / eps;

            double RHS_H = 2.0*Cd - H_star*Cf2[] 
                - (2.0*H_2star + H_star*(1.0 - ibl_H[])) * (ibl_theta[]/ue_safe) * due_dxi;
            double dHdxi_star = RHS_H / (ibl_theta[] * dHstar_dH + 1e-10);
            double H_new = ibl_H[] + ibl_dt * dHdxi_star;
   
            H_tmp[] = H_new;
        }

        // 更新变量
        foreach_cache(cutcells) {
            front_tail_handling();
            ibl_theta[] = theta_tmp[];
            ibl_H[] = H_tmp[];
        }
    }

    if (i % 100 == 0) {
        FILE* f_theta = fopen("data/theta.dat", "w");
        fprintf(f_theta, "# x theta\n");
        foreach_cache(cutcells) {
            if (y > 0.)
                fprintf(f_theta, "%.6e %.6e\n", x - stag_x, ibl_theta[]);
        }
        fflush(f_theta);
        fclose(f_theta);
        system("gnuplot plot_th.gp");

        // FILE* f_delta = fopen("data/delta.dat", "w");
        // fprintf(f_delta, "# x delta\n");
        // foreach_cache(cutcells) {
        //     if (y > 0.)
        //         fprintf(f_delta, "%.6e %.6e\n", x - stag_x, delta[]);
        // }
        // fflush(f_delta);
        // fclose(f_delta);
        // system("gnuplot plot_delta.gp");

        FILE* f_Cf = fopen("data/Cf2.dat", "w");
        fprintf(f_Cf, "# x Cf2\n");
        foreach_cache(cutcells) {
            if (y > 0.)
                fprintf(f_Cf, "%.6e %.6e\n", x - stag_x, Cf2[]*2.);
        }
        fflush(f_Cf);
        fclose(f_Cf);
        system("gnuplot plot_Cf.gp");

        FILE* f_H = fopen("data/H.dat", "w");
        fprintf(f_H, "# x H\n");
        foreach_cache(cutcells) {
            if (y > 0.)
                fprintf(f_H, "%.6e %.6e\n", x - stag_x, ibl_H[]);
        }
        fflush(f_H);
        fclose(f_H);
        system("gnuplot plot_H.gp");
    }

    }

    // foreach_cache(cutcells) {
    //     front_tail_handling();

    //     delta[] = rho[] * ibl_H * ibl_theta[] * ue[];
    // }
    
    // foreach_cache(cutcells) {
    //     front_tail_handling();

    //     coord n, b;
    //     double area = embed_geometry(point, &b, &n);
    //     foreach_dimension()
    //         n.x *= -1.;

    //     double x0 = x, y0 = y;
    //     double delta0 = delta[];
    //     double ex0 = ue_vec.x[], ey0 = ue_vec.y[];

    //     double sum_ddelta = 0., w_ddelta = 0.;
    //     foreach_neighbor(1) {
    //         if (!is_cutcell) continue;
    //         double dx = x - x0, dy = y - y0;
    //         double dist2 = sq(dx) + sq(dy);
    //         if (dist2 < 1e-10) continue;
    //         double dot = dx * ex0 + dy * ey0;
    //         if (dot < 0.) {
    //             double ds = -dot;
    //             double wt = ds / dist2;
    //             double slope = (delta0 - delta[]) / ds;
    //             sum_ddelta += wt * slope;
    //             w_ddelta += wt;
    //         }
    //     }
    //     double ddelta_ds = (w_ddelta > 0.) ? sum_ddelta / w_ddelta : 0.;
    //     double vn = ddelta_ds / rho[];
    //     double ue_local = fmax(ue[], 1e-8);
    //     double vn_max = 0.1 * ue_local;
    //     vn = fmax(fmin(vn, vn_max), -vn_max);
        
    //     foreach_dimension()
    //         w.x[] += n.x * vn * dt * area / (cm[] * Delta);
    // }
}

event output(t += 0.01) {
    FILE* f_ue = fopen("data/ue.dat", "w");
    fprintf(f_ue, "# x ue\n");
    foreach_cache(cutcells) {
        if (y > 0.)
            fprintf(f_ue, "%.6e %.6e\n", x - stag_x, ue[]/0.6);
    }
    fflush(f_ue);
    fclose(f_ue);
    system("gnuplot plot_ue.gp");
}

event end(t = 10.0);