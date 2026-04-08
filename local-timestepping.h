#include "utils.h"

extern scalar * evolving;
double (* update) (scalar * evolving, scalar * updates, double dtmax) = NULL;
double (* gradient) (double, double, double) = minmod2;
double dt = 0.;

event defaults (i = 0) {
    for (scalar s in evolving)
        s.gradient = gradient;
}

static void advance_global (scalar * output, scalar * input, 
                            scalar * updates, double dtmax) {
    if (input != output)
        trash (output);

    foreach() {
        scalar o, i, u;
        for (o, i, u in output, input, updates)
            o[] = i[] + dt * u[];
    }
}

void (* lts_advance) (scalar * output, scalar * input,
                      scalar * updates, double dt) = advance_global;

event defaults (i = 0) {
    for (scalar s in all)
        s.gradient = gradient;
    display ("box();");
}

trace
void run() {
    t = 0., iter = 0;
    dimensional (dt == DT);
    init_grid(N);

    perf.nc = perf.tnc = 0;
    perf.gt = timer_start();

    while (events (true)) {
        scalar * updates = list_clone (evolving);
        double dtmax = update (evolving, updates, DT);
        dt = dtnext (fmax(dtmax, 1e-5));

        // 一阶 Forward Euler：用 lts_advance 推进（可以是 local dt）
        lts_advance (evolving, evolving, updates, dt);

        delete (updates);
        free (updates);
        update_perf();
        iter = inext, t = tnext;
    }
    timer_print(perf.gt, iter, perf.tnc);
    free_grid();
}