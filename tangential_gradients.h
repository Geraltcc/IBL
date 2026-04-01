void cutcell_tangential_gradient(scalar f, vector g) {
    foreach() 
        foreach_dimension()
            g.x[] = 0.;
    
    foreach_cache(cutcells) {
        coord n, b, grad = {0.};
        embed_geometry(point, &b, &n);

        foreach_dimension() {
            bool has_p = (cs[1] > 0. && cs[1] < 1.);
            bool has_m = (cs[-1] > 0. && cs[-1] < 1.);

            if (has_p && has_m) 
                grad.x = (f[1] - f[-1]) / (2.*Delta);
            else if (has_p)
                grad.x = (f[1] - f[]) / Delta;
            else if (has_m)
                grad.x = (f[] - f[-1]) / Delta;
        }

        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n.x;

        foreach_dimension()
            g.x[] = grad.x - gn * n.x;
    }
}

void cutcell_tangential_gradient_extend(scalar f, vector g) {
    scalar f_ext[];
    foreach() {
        f_ext[] = f[];
        foreach_dimension()
            g.x[] = 0.;

        if (cs[] > 0.)
            continue;

        double sum = 0.;
        int count = 0.;
        foreach_dimension() {
            if (cs[1] > 0. && cs[1] < 1.) { sum += f[1]; count++; }
            if (cs[-1] > 0. && cs[-1] < 1.) { sum += f[-1]; count++; }
        }
        if (count > 0)
            f_ext[] = sum / count;
    }

    foreach_cache(cutcells) {
        coord n, b;
        embed_geometry(point, &b, &n);

        coord grad = {0.};
        foreach_dimension()
            grad.x = center_gradient(f_ext);

        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n.x;

        foreach_dimension()
            g.x[] = grad.x - gn * n.x;
    }
    delete({f_ext});
}

void cutcell_tangential_gradient_implicit_extend(scalar f, vector g) {
    foreach()
        foreach_dimension()
            g.x[] = 0.;

    foreach_cache(cutcells) {
        coord n, b;
        embed_geometry(point, &b, &n);
        
        coord grad = {0.};

        foreach_dimension() {
            double fp = (cs[1] > 0.) ? ((cs[1] < 1.) ? f[1] : f[]) : 0.;
            double fm = (cs[-1] > 0.) ? ((cs[-1] < 1.) ? f[-1] : f[]) : 0.;
            if (fs.x[] && fs.x[1])
                grad.x = (fp - fm) / (2. * Delta);
            else if (fs.x[1])
                grad.x = (fp - f[]) / Delta;
            else if (fs.x[])
                grad.x = (f[] - fm) / Delta;
        }

        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n.x;

        foreach_dimension()
            g.x[] = grad.x - gn * n.x;
    }
}

void cutcell_tangential_gradient_lsq_forward_fixed(scalar f, vector g) {
    foreach()
        foreach_dimension()
            g.x[] = 0.;

    foreach_cache(cutcells) {
        coord n_emb, bary;
        embed_geometry(point, &bary, &n_emb);

        double x0 = x, y0 = y, f0 = f[];
#if dimension > 2
        double z0 = z;
#endif

        double a00 = 0., a01 = 0., a11 = 0.;
        double b0 = 0., b1 = 0.;
#if dimension > 2
        double a02 = 0., a12 = 0., a22 = 0.;
        double b2 = 0.;
#endif
        int nn = 0;

        foreach_neighbor(1) {
            if (is_cutcell) {
                double dx = x - x0, dy = y - y0;
#if dimension > 2
                double dz = z - z0;
                double dist2 = dx*dx + dy*dy + dz*dz;
#else
                double dist2 = dx*dx + dy*dy;
#endif
                if (dist2 < 1e-20) continue;

                double w = 1. / dist2;  // 反距离平方权重
                double df = f[] - f0;

                a00 += w*dx*dx;  a01 += w*dx*dy;  a11 += w*dy*dy;
                b0  += w*dx*df;  b1  += w*dy*df;
#if dimension > 2
                a02 += w*dx*dz;  a12 += w*dy*dz;  a22 += w*dz*dz;
                b2  += w*dz*df;
#endif
                nn++;
            }
        }

        coord grad = {0.};

        // 自适应正则化：基于法向量方向
        double eps_base = 0.1 * (a00 + a11) / (nn > 0 ? nn : 1);
        
#if dimension == 2
        if (nn >= 2) {
            // 根据法向量分量调整正则化强度
            double eps_x = eps_base * (1.0 + 10.0 * n_emb.x * n_emb.x);
            double eps_y = eps_base * (1.0 + 10.0 * n_emb.y * n_emb.y);
            
            a00 += eps_x;
            a11 += eps_y;
            
            double det = a00*a11 - a01*a01;
            if (fabs(det) > 1e-30) {
                grad.x = ( a11*b0 - a01*b1) / det;
                grad.y = (-a01*b0 + a00*b1) / det;
            }
        }
#elif dimension == 3
        if (nn >= 3) {
            // 根据法向量分量调整各方向正则化
            double eps_x = eps_base * (1.0 + 10.0 * n_emb.x * n_emb.x);
            double eps_y = eps_base * (1.0 + 10.0 * n_emb.y * n_emb.y);
            double eps_z = eps_base * (1.0 + 10.0 * n_emb.z * n_emb.z);
            
            a00 += eps_x;
            a11 += eps_y;
            a22 += eps_z;

            double c00 = a11*a22 - a12*a12;
            double c01 = a02*a12 - a01*a22;
            double c02 = a01*a12 - a02*a11;
            double c11 = a00*a22 - a02*a02;
            double c12 = a01*a02 - a00*a12;
            double c22 = a00*a11 - a01*a01;
            double det = a00*c00 + a01*c01 + a02*c02;
            
            if (fabs(det) > 1e-30) {
                grad.x = (c00*b0 + c01*b1 + c02*b2) / det;
                grad.y = (c01*b0 + c11*b1 + c12*b2) / det;
                grad.z = (c02*b0 + c12*b1 + c22*b2) / det;
            }
        }
#endif

        // 投影到切平面
        double gn = 0.;
        foreach_dimension()
            gn += grad.x * n_emb.x;
        foreach_dimension()
            g.x[] = grad.x - gn * n_emb.x;
    }
}

void cutcell_tangential_gradient_lsq_backward(scalar f, vector g) {
    foreach()
        foreach_dimension()
            g.x[] = 0.;

    foreach_cache(cutcells) {
        coord n_emb, bary;
        embed_geometry(point, &bary, &n_emb);

        // ---- 构造切平面正交基 t1, t2 ----
        coord t1;
    #if dimension > 2
        coord t2;
    #endif
    
#if dimension == 2
        t1 = (coord){-n_emb.y, n_emb.x};
#elif dimension == 3
        // 选择与 n 最不平行的坐标轴做叉积
        coord ref;
        if (fabs(n_emb.x) <= fabs(n_emb.y) && fabs(n_emb.x) <= fabs(n_emb.z))
            ref = (coord){1., 0., 0.};
        else if (fabs(n_emb.y) <= fabs(n_emb.z))
            ref = (coord){0., 1., 0.};
        else
            ref = (coord){0., 0., 1.};
        // t1 = ref × n (再归一化)
        t1.x = ref.y*n_emb.z - ref.z*n_emb.y;
        t1.y = ref.z*n_emb.x - ref.x*n_emb.z;
        t1.z = ref.x*n_emb.y - ref.y*n_emb.x;
        double len1 = sqrt(t1.x*t1.x + t1.y*t1.y + t1.z*t1.z);
        t1.x /= len1;  t1.y /= len1;  t1.z /= len1;
        // t2 = n × t1
        t2.x = n_emb.y*t1.z - n_emb.z*t1.y;
        t2.y = n_emb.z*t1.x - n_emb.x*t1.z;
        t2.z = n_emb.x*t1.y - n_emb.y*t1.x;
#endif

        double x0 = x, y0 = y, f0 = f[];
#if dimension > 2
        double z0 = z;
#endif

        // ---- 2×2 正规方程 (2D 中退化为 1×1) ----
        double m11 = 0.;
        double rhs1 = 0.;


    #if dimension > 2
        double m12 = 0., m22 = 0.;
        double rhs2 = 0.;
    #endif

        foreach_neighbor(1) {

            if (is_cutcell) {

                double dx = x - x0, dy = y - y0;
    #if dimension > 2
                double dz = z - z0;
                double dist2 = dx*dx + dy*dy + dz*dz;
    #else
                double dist2 = dx*dx + dy*dy;
    #endif
                if (dist2 < 1e-20) continue;

                double w = 1. / dist2;
                double df = f[] - f0;

                // 投影到切平面基向量
                double s1 = dx*t1.x + dy*t1.y;

    #if dimension > 2
                s1 += dz*t1.z;

                double s2 = dx*t2.x + dy*t2.y;
                s2 += dz*t2.z;

                m12 += w*s1*s2;  m22 += w*s2*s2;
                rhs2 += w*df*s2;
    #endif
                m11 += w*s1*s1;  
                rhs1 += w*df*s1;  
            }
        }

        // ---- 求解 ----
        double g1 = 0.;
    #if dimension > 2
        double g2 = 0.;
    #endif

#if dimension == 2
        if (fabs(m11) > 1e-30)
            g1 = rhs1 / m11;
#elif dimension == 3

        double det = m11*m22 - m12*m12;
        if (fabs(det) > 1e-30) {
            g1 = ( m22*rhs1 - m12*rhs2) / det;
            g2 = (-m12*rhs1 + m11*rhs2) / det;
        }
       
#endif

        // ---- 反投影到笛卡尔坐标 ----
        g.x[] = g1*t1.x;
#if dimension > 1
        g.y[] = g1*t1.y;
#endif
#if dimension > 2
        g.x[] += g2 * t2.x;
        g.y[] += g2 * t2.y;
        g.z[] = g1*t1.z + g2*t2.z;
#endif
    }
}