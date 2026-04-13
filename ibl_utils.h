double cal_Hk(double Ma, double H) {
    double Hk = (H - 0.290 * sq(Ma)) / (1. + 0.113 * sq(Ma));
    return fmax(Hk, 0.);
}

double cal_Cf2(double Ma, double Re, double Hk) {
    double Fc = sqrt(1. + 0.2*sq(Ma));
    return 0.5 * fmax(0.3 * exp(-1.33 * Hk) 
        * pow(log10(fmax(Re/Fc, 2.0)), -1.74 - 0.31*Hk)
        + 0.00011 * (tanh(4. - Hk/0.875) - 1.), 0.);
}

/*
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
*/

double cal_H_star(double Re, double Hk) {
    double H0 = (Re > 400.) ? (3. + 400. / Re) : 4.;
    double Re_limit = fmax(Re, 200.);
    double base = 1.505 + 4./Re_limit;
    double term = 0.;
    if (Hk < H0) {
        term = (2. - 1.505 - 4./Re_limit) * sq((H0 - Hk) / (H0 - 1.)) * 1.5/(Hk + 0.5);
    } else {
        double dH = Hk - H0;
        double log_Re = log10(Re);
        if (log_Re < 0.1) log_Re = 0.1;
        term = sq(dH) * (0.04/Hk + 0.007 * log_Re/sq(dH + 4./log_Re));
    }
    return base + term;
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

double cal_delta(double theta, double Hk, double delta_star) {
    double delta = theta * (3.15 + 1.72 / (Hk - 1.)) + delta_star;
    // double delta = theta * (3.15 + 1.72 / (Hk - 1.));
    return delta;
}

double cal_Ctau_eq(double H, double Hk, double H_star, double Us) {
    double Ctau_eq = H_star * 0.015 / (1. - Us) * cube(Hk - 1.) / sq(Hk) / H;
    return fmax(Ctau_eq, 0.); 
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

double upwind_value(scalar s, coord e, Point point) {
    coord x0 = {x, y};
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
            sum += wt * s[];
            w += wt;
        }
    }
    return (w > 0.) ? sum / w : 0.;
}

double upwind_arc(coord e, Point point) {
    coord x0 = {x, y};
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
            sum += wt * sqrt(dist2);
            w += wt;
        }
    }
    return (w > 0.) ? sum / w : 0.;
}