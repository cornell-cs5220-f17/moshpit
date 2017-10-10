# Collective motion at heavy metal concerts

This project was adapted from a code by Matthey Bierbaum, which
implements a simple agent model of collective motion in mosh pits.
Your goal is to understand the performance of the code, tune and
parallelize it, and play with it in some interesting way.  More
specifically, you should do the following:

- Tune the serial code and report on any bottlenecks.  Is the time
  going where you expected?  Ideally, both your serial and parallel
  performance tuning should be informed by profiling experiments.
- Parallelize your code using OpenMP, and spend some time tuning the
  parallel implementation.  You will need to go beyond just adding
  pragmas before your `for` loops to get reasonable speed!
- Profile the code and set up strong and weak scaling studies.  What is the
  approximate cost per particle?  Per frame?  Does it scale linearly?
  How large a system is tractable?  What is the speedup?
- If you have extra time, play a little!  Improve or extend the code
  in some way that appeals to you.  If this project appeals to you,
  extensions could easily be the basis of a class project.

*Note*: The default visualizer uses a CSV text output; if you would
like, you may switch to the binary output format that you should be
able to view using [this Javascript viewer](http://www.cs.cornell.edu/courses/cs5220/2017fa/viz.html).

## Resources

- The [annotated code](code.md)
- [Matt's original code](https://github.com/mattbierbaum/moshpits)
- A [JavaScript version!](http://github.com/mattbierbaum/moshpits.js)
- [The collective motion project page](http://cohengroup.lassp.cornell.edu/projects/collective-motion-mosh-pits)
- My [S14 particle project](https://bitbucket.org/dbindel/cs5220-s14/wiki/HW3)
- The [OpenMP tutorials](http://www.openmp.org/resources/tutorials-articles/)
