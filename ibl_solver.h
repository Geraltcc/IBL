#include "ibl_utils.h"

#define front_tail_handling() \
    double dist_stag = sqrt(sq(x - stag_x) + sq(y - stag_y)); \
    if (dist_stag < 3.0 * Delta) { \
        ibl_theta[] = 1e-6; \
        ibl_H[] = 1.5; \
        ibl_Ctau[] = 0.005; \
        continue; \
    } \

#define MAX_LOCAL_ITER 20
#define LOCAL_TOL 1e-10

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

        double slope_limit = 0.5 * ue_safe / dxi;
        due_dxi = fmax(fmin(due_dxi, slope_limit), -slope_limit);


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

    double dtau = 1.;
    foreach_cache(cutcells) {
        front_tail_handling();

        coord e = {ue_vec.x[], ue_vec.y[]};
        
        double theta = ibl_theta[];
        double H = ibl_H[];
        double H_old = H;

        double due_dxi = central_grad(ue, e, point);
        double ue_safe = fmax(ue[], 1e-8);
        double dxi = upwind_arc(e, point);

        double slope_limit = 0.5 * ue_safe / dxi;
        due_dxi = fmax(fmin(due_dxi, slope_limit), -slope_limit);

        double H_up = upwind_value(ibl_H, e, point);
        double Re_th = fmax(rho[] * ue_safe*theta/ibl_nu, 10.);

        
        for (int m = 0; m < MAX_LOCAL_ITER; m++) {
            double Hk = fmax(cal_Hk(Me[], H), 1.05);
            double Hstar = cal_H_star(Re_th, Hk);
            double H2star = cal_H_2star(Me[], Hk);
            double Cf2 = cal_Cf2(Me[], Re_th, Hk);
            double Us = cal_Us(H, Hk, Hstar);
            double Ctau = ibl_Ctau[];
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
            double Cd_pp = cal_Cd(Cf2_pp, Us_pp, Ctau);
           
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
            H = fmin(fmax(H, 1.05), 5.0);
        }
        ibl_H[] = H;
    }
    boundary({ibl_H});

    foreach_cache(cutcells) {
        delta_star[] = ibl_H[] * ibl_theta[];
    }

    double dtau_C = 1e-2;
    foreach_cache(cutcells) {
        front_tail_handling();

        coord e = {ue_vec.x[], ue_vec.y[]};
        double ue_safe = fmax(ue[], 1e-8);

        double theta = ibl_theta[];
        double H = ibl_H[];
        double Ctau = ibl_Ctau[];
        double Ctau_old = Ctau;

        double Ctau_up = upwind_value(ibl_Ctau, e, point);
        double dxi = upwind_arc(e, point);

        double Re_th = fmax(rho[]*ue_safe*theta/ibl_nu, 10.);
        double Hk = fmax(cal_Hk(Me[], H), 1.05);
        double Hstar = cal_H_star(Re_th, Hk);
        double Us = cal_Us(H, Hk, Hstar);
        double Ctau_eq = cal_Ctau_eq(H, Hk, Hstar, Us);
        double delta = cal_delta(theta, Hk, delta_star[]);

        for (int m = 0; m < MAX_LOCAL_ITER; m++) {
            double Ctau_safe = fmax(Ctau, 1e-12);
            double delta_safe = fmax(delta, 1e-12);

            double S_Ctau = 4.2*Ctau_safe/delta_safe
                           * (sqrt(fmax(Ctau_eq, 0.)) - sqrt(Ctau_safe));

            double r_C = (Ctau - Ctau_old) / dtau_C
                       + (Ctau - Ctau_up) / dxi
                       - S_Ctau;

            if (fabs(r_C) < LOCAL_TOL) break;

            double eps_C = fmax(1e-6*fabs(Ctau), 1e-10);
            double Ctau_p = Ctau + eps_C;
            double S_Ctau_p = 4.2 * Ctau_p / delta_safe
                            * (sqrt(fmax(Ctau_eq, 0.)) - sqrt(Ctau_p));

            double r_C_p = (Ctau_p - Ctau_old) / dtau_C
                       + (Ctau_p - Ctau_up) / dxi
                       - S_Ctau_p;

            double dr_dCtau = (r_C_p - r_C) / eps_C;
            double delta_Ctau = -r_C / dr_dCtau;

            if (delta_Ctau > 0.3) delta_Ctau = 0.3;
            if (delta_Ctau < -0.3) delta_Ctau = -0.3;

            Ctau += delta_Ctau;
            Ctau = fmin(fmax(Ctau, 1e-12), 0.3);
        }
        ibl_Ctau[] = Ctau;
    }
    boundary({ibl_Ctau});
}