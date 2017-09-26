#ifndef PHYSICS_H
#define PHYSICS_H

//ldoc
/**
 * # Physics interfaces
 *
 * There are a variety of parameters that characterize the physics,
 * and it is worthwhile to keep them all in one place.  We use the
 * `set_default_params` call to set default values for each of these
 * parameters; those defaults may then be over-ridden.
 */

typedef struct sim_param_t {
    float radius;  // Interaction radius for contact forces
    float epsilon; // Contact force multiplier
    float vhappy;  // Preferred speed for active moshers
    float alpha;   // Flocking force multiplier
    float sigma;   // Noise variance
    float damp;    // Damping force multiplier
} sim_param_t;

void set_default_params(sim_param_t* p);

/**
 * The main steps in the "physics" part of the code are computing
 * particle interaction forces (`compute_forces`) and advancing the
 * particles in time with a "leapfrog" step (aka Newton-Stormer-Verlet).
 */

typedef struct particles_t particles_t;
void compute_forces(particles_t* particles, sim_param_t* params);
void leapfrog_step(particles_t* restrict particles, float dt);

//ldoc off
#endif /* PHYSICS_H */
