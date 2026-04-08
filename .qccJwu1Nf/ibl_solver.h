#ifndef BASILISK_HEADER_1
#define BASILISK_HEADER_1
#line 1 "./ibl_solver.h"
#include "ibl_utils.h"

#define front_tail_handling() \
    double dist_stag = sqrt(sq(x - stag_x) + sq(y - stag_y)); \
    double dist_tail = fabs(x - tail_x); \
    if (dist_stag < 3.0 * Delta || dist_tail < 0.01) { \
        ibl_theta[] = 1e-6; \
        ibl_H[] = 2.2216; \
        continue; \
    } \

#define MAX_LOCAL_ITER 20
#define LOCAL_TOL 1e-10

#define ibl_nu 1e-7
#define ibl_H_temp 1.8

scalar Re_th[];
scalar Hk[];
scalar H_star[];

void ibl_solver() {

    foreach_cache(cutcells) {
        Re_th[] = fmax(rho[] * fmax(ue[], 1e-8)*ibl_theta[]/ibl_nu, 10.);
        Hk[] = fmax(cal_Hk(Me[], ibl_H[]), 1.05);
        H_star[] = cal_H_star(Re_th[], Hk[]);
    }
    boundary({Re_th, Hk, H_star});

    foreach_cache(cutcells) {
        front_tail_handling();

        coord e = {ue_vec.x[], ue_vec.y[]};
        double due_dxi = central_grad(ue, e, point);
        double ue_safe = fmax(ue[], 1e-8);

        double theta_up = upwind_value(ibl_theta, e, point);
        double dxi = upwind_arc(e, point);

        // local Newton sigma*Q - S(Q) = b
        double theta = ibl_theta[];
        double H = ibl_H[];
        double Cf2_val = 0.;

        for (int m = 0; m < MAX_LOCAL_ITER; m++) {
            double Re_th_cur = fmax(rho[] * ue_safe * theta / ibl_nu, 1e-4);
            Cf2_val = cal_Cf2(Me[], Re_th_cur, Hk[]);

            double S_theta = Cf2_val - (2.+H-sq(Me[]))*theta/ue_safe*due_dxi;
            double r_theta = (theta - theta_up) / dxi - S_theta;

            if (fabs(r_theta) < LOCAL_TOL) break;

            double eps = fmax(1e-6 * fabs(theta), 1e-10);
            double theta_p = theta + eps;
            double Re_th_p = fmax(rho[]*ue_safe*theta_p/ibl_nu, 1e-4);
            double Cf2_val_p = cal_Cf2(Me[], Re_th_p, Hk[]);
            double S_theta_p = Cf2_val_p - (2.+H-sq(Me[]))*theta_p/ue_safe*due_dxi;
            double r_theta_p = (theta_p - theta_up) / dxi - S_theta_p;

            double dr_dtheta = (r_theta_p - r_theta) / eps;
            double delta_theta = -r_theta / dr_dtheta;

            theta += delta_theta;
        }
        Cf2[] = Cf2_val;
        ibl_theta[] = theta;
    }
    boundary({ibl_theta});

    // ===================== H equation ===========================

    double dtau = 1e-2;
    foreach_cache(cutcells) {
        front_tail_handling();

        coord e = {ue_vec.x[], ue_vec.y[]};
        
        double theta = ibl_theta[];
        double H = ibl_H[];
        double H_old = H;

        double due_dxi = central_grad(ue, e, point);
        double ue_safe = fmax(ue[], 1e-8);
        double H_up = upwind_value(ibl_H, e, point);
        double dxi = upwind_arc(e, point);
        double Re_th = fmax(rho[] * ue_safe*theta/ibl_nu, 10.);

        
        for (int m = 0; m < MAX_LOCAL_ITER; m++) {
            double Hk = fmax(cal_Hk(Me[], H), 1.05);
            double Hstar = cal_H_star(Re_th, Hk);
            double H2star = cal_H_2star(Me[], Hk);
            double Cf2 = cal_Cf2(Me[], Re_th, Hk);
            double Us = cal_Us(H, Hk, Hstar);
            double Ctau = cal_Ctau_eq(H, Hk, Hstar, Us);
            double Cd = cal_Cd(Cf2, Us, Ctau);

            double RHS = 2.*Cd - Hstar*Cf2
                - (2.*H2star + Hstar*(1. - H)) * theta/ue_safe*due_dxi;

            // dH*/dH
            double eps_h = fmax(1e-6*fabs(H), 1e-10);
            double Hk_p = fmax(cal_Hk(Me[], H + eps_h), 1.05);
            double Hstar_p = cal_H_star(Re_th, Hk_p);
            double dHs_dH = (Hstar_p - Hstar) / eps_h;

            double denom = theta * dHs_dH;
            if (fabs(denom) < 1e-15) denom = (denom < 0.) ? -1e-15 : 1e-15;
            double S_H = RHS / denom;

            double r_H = (H - H_old) / dtau
                       + (H - H_up) / dxi
                       - S_H;

            if (fabs(r_H) < LOCAL_TOL) break;

            double H_pp = H + eps_h;
            double Hk_pp = fmax(cal_Hk(Me[], H_pp), 1.05);
            double Hstar_pp = cal_H_star(Re_th, Hk_pp);
            double H2star_pp = cal_H_2star(Me[], Hk_pp);
            double Cf2_pp = cal_Cf2(Me[], Re_th, Hk_pp);
            double Us_pp = cal_Us(H_pp, Hk_pp, Hstar_pp);
            double Ctau_pp = cal_Ctau_eq(H_pp, Hk_pp, Hstar_pp, Us_pp);
            double Cd_pp = cal_Cd(Cf2_pp, Us_pp, Ctau_pp);
           
            double RHS_pp = 2.*Cd_pp - Hstar_pp*Cf2_pp
                - (2.*H2star_pp + Hstar_pp*(1. - H_pp)) * theta/ue_safe*due_dxi;

            double dHsdH_pp = (cal_H_star(Re_th, fmax(cal_Hk(Me[],H_pp+eps_h),1.05)) 
                - Hstar_pp) / eps_h;

            double denorm_pp = theta * dHsdH_pp;
            if (fabs(denorm_pp) < 1e-15) denorm_pp = (denorm_pp < 0.) ? -1e-15 : 1e-15;
            double S_H_pp = RHS_pp / denorm_pp;

            double r_H_pp = (H_pp - H_old) / dtau
                       + (H_pp - H_up) / dxi
                       - S_H_pp;

            double dr_dH = (r_H_pp - r_H) / eps_h;
            double delta_H = -r_H / dr_dH;

            if (delta_H > 0.3) delta_H = 0.3;
            if (delta_H < -0.3) delta_H = -0.3;

            H += delta_H;
            H = fmin(fmax(H, 1.05), 3.0);
        }
        ibl_H[] = H;
    }
    boundary({ibl_H});

    foreach_cache(cutcells) {
        delta[] = ibl_H[] * ibl_theta[];
    }
}
#endif
