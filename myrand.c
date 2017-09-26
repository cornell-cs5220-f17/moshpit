#include "myrand.h"

static unsigned long long int vseed;
static unsigned long long int vran;


void ran_seed(long j)
{
    vseed = j;  vran = 4101842887655102017LL;
    vran ^= vseed;
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    vran = vran * 2685821657736338717LL;
}


float ran_ran2()
{
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    unsigned long long int t = vran * 2685821657736338717LL;
    return 5.42101086242752217e-20*t;
}
