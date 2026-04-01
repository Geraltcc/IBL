typedef struct {
    int i, j, level;
    double x, y;
    double s;
    double ue, Me, theta, Cf2;
} SurfaceNode;

#define MAX_SURFACE_NODES 2048
SurfaceNode chain_upper[MAX_SURFACE_NODES];
SurfaceNode chain_lower[MAX_SURFACE_NODES];
int n_upper = 0, n_lower = 0;

scalar visited[];

void build_surface_chains() {
    foreach()
        visited[] = 0.;

    Point stag_point = {0.};
    double ue_min_local = HUGE;
    foreach_cache(cutcells) {
        if (ue[] < ue_min_local) {
            ue_min_local = ue[];
            stag_point = point;
        }
    }

    coord n_stag, b_stag;
    embed_geometry(stag_point, &b_stag, &n_stag);

    coord t_pos = { -n_stag.y,  n_stag.x };
    coord t_neg = {  n_stag.y, -n_stag.x };
}