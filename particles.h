#ifndef PARTICLES_H
#define PARTICLES_H

#include <stdio.h>

//ldoc on
/**
 * # Particle interfaces
 *
 * Agents (particles) in the simulation are characterized by a
 * position, velocity, and type: "black" moshers are passive and
 * prefer to stay in one place; "red" moshers want to move around.
 * Positions are in an $L$-by-$L$ periodic domain, which is split into
 * bins (`nbinx` on a side) such that we only need to check for agent
 * interactions across neighboring bins.
 *
 * The particle positions and velocities (and accumulated forces acting
 * on the particles) are kept in parallel arrays, each consisting of
 * $x$ and $y$ components for each particle in turn.  Another array
 * accounts for the particle types.  The assignment of particles to
 * bins is done via a linked list structure: the `nbinx`-by-`nbinx`
 * array `cells` gives the index of the particle at the head of the list
 * for each of the cells, and the `next` array for each particles indicates
 * whatever is next in the list in the same cell.  A value of `-1` indicates
 * the end of the list.
 */

#define BLACK   0    // Passive mosher tag
#define RED     1    // Active mosher tag

typedef struct particles_t {
    int nbinx;           // Number of bins per direction
    int N;               // Number of particles
    float L;             // Domain size
    float* restrict x;   // Positions
    float* restrict v;   // Velocities
    float* restrict f;   // Forces
    int* restrict type;  // Particle types
    int* restrict cells; // Neighbor lists
    int* restrict next;  // Next pointers
} particles_t;

/**
 * ## Particle structure memory management
 */

particles_t* alloc_particles_t(int nbinx, int N, float L);
void         free_particles_t(particles_t* particles);

/**
 * ## Particle structure initialization
 *
 * When we set up a simulation, we need to give initial conditions for
 * the particle locations and velocities.  The `init_circle` function
 * initializes the particles into a rotating circle initially; the
 * `init_ric` gives an alternate initialization.  In both cases, the 
 * `speed` parameter indicates the preferred speed of the "red" (active)
 * moshers.
 */

void init_ric   (particles_t* particles, float speed);
void init_circle(particles_t* particles, float speed);

/**
 * ## Particle I/O
 *
 * Writing out the particle positions is expensive, but it seems
 * necessary if we want to make pretty pictures of how they are all
 * moving.  The particle I/O system writes out a simple text-based CSV
 * (comma-separated value) file with the fields:
 *
 * - `PTag`: Indicates whether particles is black (0) or red (1)
 * - `PId`: Unique integer identifier for the particle
 * - `PLocX`, `PLocY`: Position components
 * - `PDirX`, `PDirY`: Velocity components
 *
 * The CSV file may have several frames; they are assumed to be in
 * consecutive order.
 */

FILE* start_frames(const char* fname);
void write_frame(FILE* fp, particles_t* particles);
void end_frames(FILE* fp);

/**
 * ## Neighbor list computation
 *
 * The basic neighbor list computation is in `compute_nbr_lists`.
 * The `coord_to_index` helper function is just there to convert
 * a floating point coordinate to an integer bin index between
 * zero and `nbinx-1`.  There are other ways to assign coordinates
 * to indices -- perhaps a Z-Morton encoding would be better for
 * spatial locality, for example.
 */

void compute_nbr_lists(particles_t* particles);

inline int coord_to_index(float x, int nbinx, float L)
{
    return (int) (x/L * nbinx);
}

//ldoc off
#endif /* PARTICLES_H */
