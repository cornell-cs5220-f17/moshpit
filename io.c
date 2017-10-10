#include "particles.h"

#include <stdio.h>
#include <stdlib.h>

//ldoc
/**
 * ## I/O subsystem
 * 
 * The I/O subsystem writes out a text file with particle information;
 * given half a chance, this takes more time than just about anything else
 * in the simulation!  We could be somewhat more space efficient using
 * a binary format (and you are welcome to do so, though the visualizer
 * will need corresponding changes) -- but really, the "right" approach
 * is probably to moderate our desire to dump out too much information.
 * The pictures are pretty, but summary statistics are what the researchers
 * who designed this simulation were really after.
 */


FILE* start_frames(const char* fname, int npart)
{
    FILE* fp = fopen(fname, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s for output\n", fname);
        exit(-1);
    }
    uint16_t npart_lo = (npart & 0xFFFF);
    uint16_t npart_hi = (npart >> 16);
    fwrite(&npart_lo, sizeof(uint16_t), 1, fp);
    fwrite(&npart_hi, sizeof(uint16_t), 1, fp);
    return fp;
}


void write_frame(FILE* fp, particles_t* particles)
{
    particle_t* p = particles->p;
    int n = particles->N;
    float L = particles->L;
    for (int i = 0; i < n; ++i) {
        uint16_t x[3] = {p[i].type,
                         0xFFFF * p[i].x[0]/L,
                         0xFFFF * p[i].x[1]/L};
        fwrite(x, sizeof(uint16_t), 3, fp);
    }
}


void end_frames(FILE* fp)
{
    fclose(fp);
}

