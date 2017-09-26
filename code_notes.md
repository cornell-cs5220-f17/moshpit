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

