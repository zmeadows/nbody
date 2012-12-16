#ifndef NBODY_H
#define NBODY_H

#include "base.h"
#include "stddef.h"
#include "stdlib.h"
#include "time.h"
#include "float.h"
#include "stdio.h"
#include "math.h"

/* nbody simulator using fourth order hermite integration */
class nbody_solver {
private:
    const size_t DIM;
    const size_t N;

    const real_t dt_param;
    real_t dt;

    real_t *mass;

    real_t *pos;
    real_t *vel;
    real_t *acc;
    real_t *jerk;

    real_t *old_pos;
    real_t *old_vel;
    real_t *old_acc;
    real_t *old_jerk;

    real_t time_elapsed;
    size_t steps_taken;
    real_t potential;
    real_t collision_time;

    void predict_step(void);
    void advance_step(void);
    void correct_step(void);


public:
    nbody_solver(
            const size_t _DIM,
            const size_t _N,
            const real_t _dt_param,
            real_t *init_pos,
            real_t *init_vel,
            real_t *init_mass
            );

    nbody_solver(
            const size_t _DIM,
            const size_t _N,
            const real_t _dt_param
            );
    ~nbody_solver(void);
    
    real_t get_time_elapsed(void) const { return time_elapsed; };
    size_t get_steps_taken(void) const { return steps_taken; };
    size_t get_num_bodies(void) const { return N; };
    size_t get_dimension(void) const { return DIM; };

    void evolve_step(void);
    
    void print_stats(void);
    void print_diagnostics(void);
    
    real_t get_x_pos(size_t id) const { return pos[id*DIM]; };
    real_t get_y_pos(size_t id) const { return pos[id*DIM + 1]; };
    real_t get_mass(size_t id) const { return mass[id]; };
};

#endif
