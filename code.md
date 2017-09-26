# Annotated code for the moshpit project

This project was adapted from [a code by Matthey Bierbaum][mosh-gh],
which implements a simple agent model of
[collective motion in mosh pits][mosh-page].
The [official paper] on the project was published in Physical Review
Letters; there is also a shorter [arXiv version].
Your goal is to understand the performance of the code,
tune and parallelize it, and play with it in some interesting way.

Matt's original code is already parallelized with OpenMP, and you are
welcome to use his implementation as the basis for your work (with
appropriate citation).  However, I would suggest first addressing some
of the data locality issues in the access pattern.  Also, as noted in
class, with a little tuning you are unlikely to get good speedup for
small numbers of agents if you simply use a parallel for as is done in
this code -- better speedups require some more careful thought.

[mosh-gh]: https://github.com/mattbierbaum/moshpits
[mosh-page]: http://cohengroup.lassp.cornell.edu/projects/collective-motion-mosh-pits
[official paper]: http://prl.aps.org/abstract/PRL/v110/i22/e228701
[arXiv version]: http://arxiv.org/abs/1302.1886

# Particle interfaces

Agents (particles) in the simulation are characterized by a
position, velocity, and type: "black" moshers are passive and
prefer to stay in one place; "red" moshers want to move around.
Positions are in an $L$-by-$L$ periodic domain, which is split into
bins (`nbinx` on a side) such that we only need to check for agent
interactions across neighboring bins.

The particle positions and velocities (and accumulated forces acting
on the particles) are kept in parallel arrays, each consisting of
$x$ and $y$ components for each particle in turn.  Another array
accounts for the particle types.  The assignment of particles to
bins is done via a linked list structure: the `nbinx`-by-`nbinx`
array `cells` gives the index of the particle at the head of the list
for each of the cells, and the `next` array for each particles indicates
whatever is next in the list in the same cell.  A value of `-1` indicates
the end of the list.
    
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
    
## Particle structure memory management
    
    particles_t* alloc_particles_t(int nbinx, int N, float L);
    void         free_particles_t(particles_t* particles);
    
## Particle structure initialization

When we set up a simulation, we need to give initial conditions for
the particle locations and velocities.  The `init_circle` function
initializes the particles into a rotating circle initially; the
`init_ric` gives an alternate initialization.  In both cases, the 
`speed` parameter indicates the preferred speed of the "red" (active)
moshers.
    
    void init_ric   (particles_t* particles, float speed);
    void init_circle(particles_t* particles, float speed);
    
## Particle I/O

Writing out the particle positions is expensive, but it seems
necessary if we want to make pretty pictures of how they are all
moving.  The particle I/O system writes out a simple text-based CSV
(comma-separated value) file with the fields:

- `PTag`: Indicates whether particles is black (0) or red (1)
- `PId`: Unique integer identifier for the particle
- `PLocX`, `PLocY`: Position components
- `PDirX`, `PDirY`: Velocity components

The CSV file may have several frames; they are assumed to be in
consecutive order.
    
    FILE* start_frames(const char* fname);
    void write_frame(FILE* fp, particles_t* particles);
    void end_frames(FILE* fp);
    
## Neighbor list computation

The basic neighbor list computation is in `compute_nbr_lists`.
The `coord_to_index` helper function is just there to convert
a floating point coordinate to an integer bin index between
zero and `nbinx-1`.  There are other ways to assign coordinates
to indices -- perhaps a Z-Morton encoding would be better for
spatial locality, for example.
    
    void compute_nbr_lists(particles_t* particles);
    
    inline int coord_to_index(float x, int nbinx, float L)
    {
        return (int) (x/L * nbinx);
    }
    

# Physics interfaces

There are a variety of parameters that characterize the physics,
and it is worthwhile to keep them all in one place.  We use the
`set_default_params` call to set default values for each of these
parameters; those defaults may then be over-ridden.
    
    typedef struct sim_param_t {
        float radius;  // Interaction radius for contact forces
        float epsilon; // Contact force multiplier
        float vhappy;  // Preferred speed for active moshers
        float alpha;   // Flocking force multiplier
        float sigma;   // Noise variance
        float damp;    // Damping force multiplier
    } sim_param_t;
    
    void set_default_params(sim_param_t* p);
    
The main steps in the "physics" part of the code are computing
particle interaction forces (`compute_forces`) and advancing the
particles in time with a "leapfrog" step (aka Newton-Stormer-Verlet).
    
    typedef struct particles_t particles_t;
    void compute_forces(particles_t* particles, sim_param_t* params);
    void leapfrog_step(particles_t* restrict particles, float dt);
    

# Random number generation

A pseudo-random number generator (PRNG) is a function that produces
a sequence of random-looking outputs.  Typically, a PRNG consists
of a finite state machine (often implemented as some nonlinear map
involving arithmetic, bit-twiddling, and modular operations) and a
function that maps the current state to an output variate.
There is a *lot* of art that goes into pseudo-random number
generators (PRNGs), balancing cost of producing a pseudo-random
number against the quality of the number.  For Monte Carlo
computations, quality is measured in terms of various statistical
tests -- the bar (and stakes) are higher when pseudo-random number
generators are used for cryptographic applications.

Parallel pseudo-random number generation is tricky in that
we don't want the streams to be correlated.  We also typically want
to use a generator which avoids global state (and is thus thread
safe).  The generator code here, inherited from the original mosh
pit code, appears to be the Numerical Recipes generator Ranq1.
As implemented, it is not even thread-safe!  You may decide you 
don't care -- but if you decide to care, you are advised to probably
use one of the MKL parallel PRNGs.

    void  ran_seed(long j);
    float ran_ran2();
    

# Particle data stucture implementation

The particle data structure is one of the places where you may
want to improve the implementation.  There are several aspects
that are not entirely satisfactory; are we better off with the
parallel array layout we have now (which is good for vectorization,
but perhaps not locality), or with an array-of-structs layout?
Should the particle type (which is really only a single bit of
information) be stored as a full integer, or as some smaller data
type?  Is there a better way of handling the bin data structures
(noting that this is not the bin data structure from the original
code)?  Play with things and find out!

## Allocate and free particle data
    
    particles_t* alloc_particles_t(int nbinx, int N, float L)
    {
        int size_total = nbinx*nbinx;
        particles_t* p = malloc(sizeof(particles_t));
        p->nbinx = nbinx;
        p->N = N;
        p->L = L;
        p->type =  (int*)   calloc(  N, sizeof(int));
        p->x =     (float*) calloc(2*N, sizeof(float));
        p->v =     (float*) calloc(2*N, sizeof(float));
        p->f =     (float*) calloc(2*N, sizeof(float));
        p->next =  (int*)   calloc(  N, sizeof(int));
        p->cells = (int*)   calloc(size_total, sizeof(int));
        return p;
    }
    
    
    void free_particles_t(particles_t* p)
    {
        free(p->cells);
        free(p->next);
        free(p->f);
        free(p->v);
        free(p->type);
        free(p);
    }
    
    
## Neighbor list computations

Right now, computing a neighbor list means computing the `cells`
and `next` arrays; we don't touch the position, velocity, or force
arrays.  But perhaps we should!  As currently written, traversing
a neighbor list involves jumping all over memory to grab the information
for different particles; you may want to consider something that makes
better use of cache locality.
    
    void compute_nbr_lists(particles_t* particles)
    {
        int* restrict cells = particles->cells;
        int* restrict next  = particles->next;
        float* restrict x = particles->x;
        float L = particles->L;
        int N = particles->N;
        int nbinx = particles->nbinx;
    
        // Recompute neighbor list
        const int size_total = nbinx*nbinx;
        memset(cells, -1, size_total * sizeof(int));
        for (int i=0; i<N; i++) {
            int binx = coord_to_index(x[2*i+0], nbinx, L);
            int biny = coord_to_index(x[2*i+1], nbinx, L);
            int t = binx + biny*nbinx;
            next[i] = cells[t];
            cells[t] = i;
        }
    }
    
    
## Initialization

The `initial_circle` function puts all the active moshers in a big
circle in the middle of the domain at the start; the `init_ric`
function distributes the active moshers uniformly through the
domain and assigns them a random velocity as well as position.
It might be interesting to see how the behavior depends on this
initial distribution.
    
    
    void init_ric(particles_t* particles, float speed)
    {
        float* restrict x = particles->x;
        float* restrict v = particles->v;
        int* restrict type = particles->type;
        int N = particles->N;
        float L = particles->L;
    
        for (int i=0; i<N; i++) {
            float t = 2*M_PI*ran_ran2();
    
            x[2*i+0] = L*ran_ran2();
            x[2*i+1] = L*ran_ran2();
    
            if (ran_ran2() > 0.16){
                v[2*i+0] = 0.0;
                v[2*i+1] = 0.0;
                type[i] = BLACK;
            } else {
                v[2*i+0] = speed * sin(t);
                v[2*i+1] = speed * cos(t);
                type[i] = RED;
            }
        }
    }
    
    
    void init_circle(particles_t* particles, float speed)
    {
        float* restrict x = particles->x;
        float* restrict v = particles->v;
        int* restrict type = particles->type;
        int N = particles->N;
        float L = particles->L;
        
        for (int i=0; i<N; i++){
            float tx = L*ran_ran2();
            float ty = L*ran_ran2();
            float tt = 2*M_PI*ran_ran2();
    
            x[2*i+0] = tx;
            x[2*i+1] = ty;
    
            // the radius for which 30% of the particles are red on avg
            float dd2 = (tx-L/2)*(tx-L/2) + (ty-L/2)*(ty-L/2);
            float rad2 = 0.16*L*L / M_PI;
    
            if (dd2 < rad2)
                type[i] = RED;
            else
                type[i] = BLACK;
    
            if (type[i] == RED) {
                v[2*i+0] = speed*cos(tt);
                v[2*i+1] = speed*sin(tt);
            } else {
                v[2*i+0] = 0.0;
                v[2*i+1] = 0.0;
            }
        }
    }
    
    
## I/O subsystem

The I/O subsystem writes out a text file with particle information;
given half a chance, this takes more time than just about anything else
in the simulation!  We could be somewhat more space efficient using
a binary format (and you are welcome to do so, though the visualizer
will need corresponding changes) -- but really, the "right" approach
is probably to moderate our desire to dump out too much information.
The pictures are pretty, but summary statistics are what the researchers
who designed this simulation were really after.
    
    
    FILE* start_frames(const char* fname)
    {
        FILE* fp = fopen(fname, "w");
        if (fp == NULL) {
            fprintf(stderr, "Could not open %s for output\n", fname);
            exit(-1);
        }
        fprintf(fp, "PTag,PId,PLocX,PLocY,PDirX,PDirY\n");
        return fp;
    }
    
    
    void write_frame(FILE* fp, particles_t* particles)
    {
        float* x = particles->x;
        float* v = particles->v;
        int* tag = particles->type;
        int n = particles->N;
        float L = particles->L;
        for (int i = 0; i < n; ++i) {
            fprintf(fp, "%d,%d,%g,%g,%g,%g\n", tag[i], i+1,
                    x[2*i+0]/L, x[2*i+1]/L,
                    v[2*i+0]/L, v[2*i+1]/L);
        }
    }
    
    
    void end_frames(FILE* fp)
    {
        fclose(fp);
    }
    


# Main driver

The simulation takes a list of physical parameters, an initial
condition, a number of particles, and some information about
the time stepping and the rate of frame output for visualization.
    
    void simulate(sim_param_t* params,   // Physics parameters
                  const char* ic_name,   // IC name (ric or circle)
                  int N,                 // Number of particles
                  int rendert,           // Time steps between frame dumps
                  float tfinal,          // End time
                  float dt)              // Time step
    {
        float radius     = params->radius;
        float L          = 1.03*sqrtf(M_PI*radius*radius*N);
        int   nbinx      = (int) (L/(4*radius));
    
        particles_t* particles = alloc_particles_t(nbinx, N, L);
        if (strcmp(ic_name, "ric") == 0)
            init_ric(particles, params->vhappy);
        else
            init_circle(particles, params->vhappy);
    
        FILE* fp = start_frames("particles.csv");
        int frames = 0;
        for (float t=0.0; t < tfinal; t += dt) {
            if (rendert > 0 && frames % rendert == 0)
                write_frame(fp, particles);
            compute_nbr_lists(particles);
            compute_forces(particles, params);
            leapfrog_step(particles, dt);
            ++frames;
        }
        end_frames(fp);
        free_particles_t(particles);
    }
    
Our `main` function uses the `getopt` library to parse options,
then runs the simulation.  This is the type of thing that could
easily be done by access to a little language for reading configuration
files as well -- I like Lua for this purpose.
    
    int main(int argc, char **argv)
    {
        sim_param_t params;
        set_default_params(&params);
        const char* ic_name = "circle";
        
        int seed     = 0;      // Random seed
        int N        = 1000;   // Number of particles
        int rendert  = 0;      // Time between dumped frames
        float tfinal = 1e3;    // Final time
        float dt = 1e-1;       // Time step
    
        // Argument processing
        int c;
        extern char* optarg;
        while ((c = getopt(argc, argv, "?r:e:v:a:s:d:N:T:h:f:")) != -1) {
            switch (c) {
            case '?':
                fprintf(stderr,
                        "%s\n"
                        "\t-?: print this message\n"
                        "\t-i: initial conditions (%s)\n"
                        "\t-r: interaction radius (%g)\n"
                        "\t-e: interaction strength epsilon (%g)\n"
                        "\t-v: preferred speed for moshers (%g)\n"
                        "\t-a: flocking force multiplier (%g)\n"
                        "\t-s: noise variance (%g)\n"
                        "\t-d: damping force multiplier (%g)\n"
                        "\t-N: number of particles (%d)\n"
                        "\t-T: final time (%g)\n"
                        "\t-h: time step (%g)\n"
                        "\t-f: steps between output frames (%d)\n",
                        argv[0], ic_name,
                        params.radius, params.epsilon, params.vhappy,
                        params.alpha, params.sigma, params.damp,
                        N, tfinal, dt, rendert);
                return -1;
            case 'i': ic_name        = optarg;       break;
            case 'r': params.radius  = atof(optarg); break;
            case 'e': params.epsilon = atof(optarg); break;
            case 'v': params.vhappy  = atof(optarg); break;
            case 'a': params.alpha   = atof(optarg); break;
            case 's': params.sigma   = atof(optarg); break;
            case 'd': params.damp    = atof(optarg); break;
            case 'N': N              = atoi(optarg); break;
            case 'T': tfinal         = atof(optarg); break;
            case 'h': dt             = atof(optarg); break;
            case 'f': rendert        = atoi(optarg); break;
            default:
                fprintf(stderr, "Unknown option (-%c)\n", c);
                return -1;
            }
        }
    
        ran_seed(seed);
        simulate(&params, ic_name, N, rendert, tfinal, dt);
    
        return 0;
    }
    

