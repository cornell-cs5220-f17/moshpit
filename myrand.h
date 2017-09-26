#ifndef MYRAND_H
#define MYRAND_H

//ldoc on
/**
 * # Random number generation
 *
 * A pseudo-random number generator (PRNG) is a function that produces
 * a sequence of random-looking outputs.  Typically, a PRNG consists
 * of a finite state machine (often implemented as some nonlinear map
 * involving arithmetic, bit-twiddling, and modular operations) and a
 * function that maps the current state to an output variate.
 * There is a *lot* of art that goes into pseudo-random number
 * generators (PRNGs), balancing cost of producing a pseudo-random
 * number against the quality of the number.  For Monte Carlo
 * computations, quality is measured in terms of various statistical
 * tests -- the bar (and stakes) are higher when pseudo-random number
 * generators are used for cryptographic applications.
 *
 * Parallel pseudo-random number generation is tricky in that
 * we don't want the streams to be correlated.  We also typically want
 * to use a generator which avoids global state (and is thus thread
 * safe).  The generator code here, inherited from the original mosh
 * pit code, appears to be the Numerical Recipes generator Ranq1.
 * As implemented, it is not even thread-safe!  You may decide you 
 * don't care -- but if you decide to care, you are advised to probably
 * use one of the MKL parallel PRNGs.
 *
 */
void  ran_seed(long j);
float ran_ran2();

//ldoc off
#endif /* MYRAND_H */
