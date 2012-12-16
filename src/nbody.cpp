#include "nbody.h"

nbody_solver::nbody_solver(
    const size_t _DIM,
    const size_t _N,
    const real_t _dt_param,
    real_t *init_pos,
    real_t *init_vel,
    real_t *init_mass)
    : DIM(_DIM),
      N(_N),
      dt_param(_dt_param),
      mass(new real_t[_N]),
      pos(new real_t[_N*_DIM]),
      vel(new real_t[_N*_DIM]),
      acc(new real_t[_N*_DIM]),
      jerk(new real_t[_N*_DIM]),
      old_pos(new real_t[_N*_DIM]),
      old_vel(new real_t[_N*_DIM]),
      old_acc(new real_t[_N*_DIM]),
      old_jerk(new real_t[_N*_DIM]),
      time_elapsed(0),
      steps_taken(0)
{
    size_t i,k;
    for (i = 0; i < N; ++i) {
        mass[i] = init_mass[i];
        for (k = 0; k < DIM; ++k) {
            pos[i*DIM+k] = init_pos[i*DIM+k];
            vel[i*DIM+k] = init_vel[i*DIM+k];
        }
    }
    advance_step();
}

nbody_solver::nbody_solver(
    const size_t _DIM,
    const size_t _N,
    const real_t _dt_param)
    : DIM(_DIM),
      N(_N+1),
      dt_param(_dt_param),
      mass(new real_t[N]),
      pos(new real_t[N*DIM]),
      vel(new real_t[N*DIM]),
      acc(new real_t[N*_DIM]),
      jerk(new real_t[N*_DIM]),
      old_pos(new real_t[N*_DIM]),
      old_vel(new real_t[N*_DIM]),
      old_acc(new real_t[N*_DIM]),
      old_jerk(new real_t[N*_DIM]),
      time_elapsed(0),
      steps_taken(0)

{
    srand (time(NULL));

    size_t i,k;
    for (i = 0; i < _N; ++i) {
        mass[i] = rand_double(10,1000);
        for (k = 0; k < DIM; ++k) {
            real_t r = rand_double(20,100);
            pos[i*DIM+k] = r;
        }
        real_t theta_v = M_PI/2 - atan(pos[i*DIM]/pos[i*DIM+1]);
        real_t mag_v = 500;
        vel[i*DIM] = signum(pos[i*DIM + 1]) * cos(theta_v) * mag_v;
        vel[i*DIM+1] = - signum(pos[i*DIM + 0]) * sin(theta_v) * mag_v;


    }

    mass[N-1] = 9999000;
    pos[DIM*(N-1)] = 60;
    pos[DIM*(N-1)+1] = 60;
    vel[DIM*(N-1)] = 0;
    vel[DIM*(N-1)+1] = 0;
    advance_step();
}

nbody_solver::~nbody_solver(void) {
    delete mass;
    delete pos;
    delete vel;
    delete acc;
    delete jerk;
    delete old_pos;
    delete old_vel;
    delete old_acc;
    delete old_jerk;
}

void nbody_solver::predict_step(void) {
    size_t i,k;
    for (i = 0;i < N; ++i) {
        for (k = 0; k < DIM; ++k) {
            pos[i*DIM+k] += vel[i*DIM+k]*dt + acc[i*DIM+k]*dt*dt/2 + jerk[i*DIM+k]*dt*dt*dt/6;
            vel[i*DIM+k] += acc[i*DIM+k]*dt + jerk[i*DIM+k]*dt*dt/2;
        }
    }
}

void nbody_solver::evolve_step(void) {
    dt = dt_param * collision_time;
    
    size_t i,k,q;
    for (i = 0; i < N; ++i) {
        for (k = 0; k < DIM; ++k) {
            q = i * DIM + k;
            old_pos[q] = pos[q];
            old_vel[q] = vel[q];
            old_acc[q] = acc[q];
            old_jerk[q] = jerk[q];
        }
    }

    predict_step();
    advance_step();
    correct_step();

    time_elapsed += dt;
    ++steps_taken;
}

void nbody_solver::advance_step(void) {
    size_t i,j,k;
    for (i = 0; i < N; ++i) {
        for (k = 0; k < DIM; ++k) {
            acc[i*DIM+k] = jerk[i*DIM+k] = 0;
        }
    }

    potential = 0;

    collision_time = DBL_MAX;
    real_t collision_estimate;

    real_t *da = new real_t[DIM];
    real_t *dj = new real_t[DIM];
    real_t *rji = new real_t[DIM];
    real_t *vji = new real_t[DIM];

    for (i = 0; i < N; ++i) {
        for (j = i + 1; j < N; ++j) {

            for (k = 0; k < DIM; ++k) {
                rji[k] = pos[j*DIM+k] - pos[i*DIM+k];
                vji[k] = vel[j*DIM+k] - vel[i*DIM+k];
            }

            real_t r2 = 0;
            real_t v2 = 0;
            real_t rv_r2 = 0;

            for (k = 0; k < DIM; ++k) {
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k];
            }

            rv_r2 /= r2;
            real_t r = sqrt(r2);
            real_t r3 = r * r2;

            potential -= mass[i] * mass[j] / r;

            for (k = 0; k < DIM; ++k) {
                da[k] = rji[k] / r3;
                dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;
            }

            for (k = 0; k < DIM; ++k) {
                acc[i*DIM+k] += mass[j] * da[k];
                acc[j*DIM+k] -= mass[i] * da[k];
                jerk[i*DIM+k] += mass[j] * dj[k];
                jerk[j*DIM+k] -= mass[i] * dj[k];
            }

            collision_estimate = (r2*r2) / (v2*v2);
            if (collision_time > collision_estimate) collision_time = collision_estimate;

            real_t da2 = 0;
            for (k = 0; k < DIM; ++k) {
                da2 += da[k] * da[k];
            }

            double mij = mass[i] + mass[j];
            da2 *= mij * mij;

            collision_estimate = r2/da2;
            if (collision_time > collision_estimate) collision_time = collision_estimate;

        }
    }

    collision_time = sqrt(sqrt(collision_time));

    delete da;
    delete dj;
    delete rji;
    delete vji;
}

void nbody_solver::correct_step(void) {
    size_t i,k,q;
    for (i = 0; i < N; ++i) {
        for (k = 0; k < DIM; ++k) {
            q = i * DIM + k;
            vel[q] = old_vel[q] + (old_acc[q] + acc[q])*dt/2 + (old_jerk[q] - jerk[q])*dt*dt/12;
            pos[q] = old_pos[q] + (old_vel[q] + vel[q])*dt/2 + (old_acc[q] - acc[q])*dt*dt/12;
        }
    }
}

void nbody_solver::print_stats(void) {
    printf("%.16g\t", time_elapsed);

    size_t i,k;
    for (i = 0; i < N; ++i) {

        printf("%.16g\t", mass[i]);

        for (k = 0; k < DIM; ++k) {
            printf("%.16g\t", pos[i*DIM+k]);
        }
        
        for (k = 0; k < DIM; ++k) {
            printf("%.16g\t", vel[k]);
        }
    }
    printf("\n");
}
