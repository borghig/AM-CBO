# AM-CBO
Numerical implementation and validation of an adaptive consensus-based algorithm for multi-objective optimization (AM-CBO).
AM-CBO makes use of a schalarization startegy to simultaneously solve several scalar parametrized sub-problems. An additional adaptive strategy aims to distribute
the particles uniformly over the image space by using energy-based measures.

Version 1.0

Date 31.07.2022

------

## R e f e r e n c e s

### [1] A consensus-based algorithm for multi-objective optimization and its mean-field description.

https://arxiv.org/abs/2203.16384

and

### [2] An adaptive consensus based method for multi-objective optimization with uniform Pareto front approximation

https://arxiv.org/abs/2208.01362

by

- Giacomo &nbsp; B o r g h i &nbsp; (RWTH Aachen University & University of Ferrara), 
- Michael &nbsp; H e r t y &nbsp; (RWTH Aachen University),
- Lorenzo &nbsp; P a r e s c h i &nbsp; (University of Ferrara)

------

## D e s c r i p t i o n

MATLAB implementation of AM-CBO for bi-objective optimization problems. The method is tested against the Lamé and DOD2K benchamark problems with arbirtary dimension of the serach space
For the reader's convenience we describe the functions and the main script in what follows:

* TestRun.m: working example of AM-CBO algorithm applied to benchmark problems considered in [2]
* problem_generator.m: generates Lamé and DOD2K benchmark problems of arbitrary search space dimension
* ReferenceSolutions.mat: pre-computed reference solutions of the test problems considered in [2]
* AM-CBO.m: algorithm implementation
* auxiliary functions 
  * compute_IGD.m: computes GD and IGD perfomance metrics given a reference solution
  * compute_potential.m: computes total potential energy of a given particles distribution, together with the associated vector field.
  * compute_yalpha.m: computed weighted averages according to the Gibbs distribution as in [2]
  * faster_pareto2.m: given a set of points, deletes the dominated ones. [credits: by Ahmed  HASSAN (ahah432@yahoo.com) and Maurel AZA-GNANDJI (my.lrichy@gmail.com)]
  * lebesgue_measure.m: computes the hypervolume contribution diversity metric [Copyright (c) 2009, Yi Cao]

  
------

## C i t a t i o n s

```bibtex
@misc{https://doi.org/10.48550/arxiv.2208.01362,
  doi = {10.48550/ARXIV.2208.01362},
  url = {https://arxiv.org/abs/2208.01362},
  author = {Borghi, Giacomo and Herty, Michael and Pareschi, Lorenzo},
  keywords = {Optimization and Control (math.OC), FOS: Mathematics, FOS: Mathematics, 35Q70, 35Q84, 35Q93, 90C29, 90C56},
  title = {An adaptive consensus based method for multi-objective optimization with uniform Pareto front approximation},
  publisher = {arXiv},
  year = {2022},
  copyright = {arXiv.org perpetual, non-exclusive license}
}

@misc{borghi2022multi,
  doi = {10.48550/ARXIV.2203.16384},
  url = {https://arxiv.org/abs/2203.16384},
  author = {Borghi, Giacomo and Herty, Michael and Pareschi, Lorenzo},
  keywords = {Optimization and Control (math.OC), FOS: Mathematics, FOS: Mathematics, 35Q70, 35Q84, 35Q93, 90C29, 90C56},
  title = {A consensus-based algorithm for multi-objective optimization and its mean-field description},
  journal = {Proceedings of, the 61st IEEE Conference on Decision and Control, to appear},
  year = {2022},
    copyright = {arXiv.org perpetual, non-exclusive license}
}
```
(Template taken from https://github.com/KonstantinRiedl/CBOGlobalConvergenceAnalysis)
