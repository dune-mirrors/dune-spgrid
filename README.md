Dune-SPGrid
===========

The DUNE module dune-spgrid provides a structured, parallel grid. Though the
dune-grid module already contains two structured grids, SGrid and YaspGrid,
both lack several important features (see below for details).


Features
--------

The following table compares the features of SPGrid to those of SGrid and
YaspGrid:

|                                       | SGrid[^1] | YaspGrid | SPGrid |
| :------------------------------------ | :-------: | :------: | :----: |
| Can communicate on codimensions       | none      | all      | all    |
| Coordinate type is template parameter | no        | yes      | yes    |
| Supports anisotropic refinement       | no        | no       | yes    |
| Supports periodic boundary conditions | no        | no[^2]   | yes    |
| Supports non-blocking communication   | no        | no       | yes    |

[^1]: SGrid is deprecated in DUNE 2.4 and will be removed in DUNE 3.0.
[^2]: YaspGrid supports a different concept of periodicity.

SPGrid supports different (global) refinement techniques, selected by a
template parameter. Some refinement techniques allow an optional parameter,
the refinement policy, to be passed to globalRefine. Currently, isotropic,
anisotropic and bisection refinement are supported. If no policy is given,
both, isotropic and anisotropic refinement, split each cube into 2^dim child
cubes. For anisotropic refinement, a policy may be used to say which directions
to split. Bisection refinement always split a cube into 2 child cubes; the
split direction can be given by the policy. If no policy is given, the split
directions are cycled through. By default, isotropic refinement is used.


Preprocessor Magic
------------------

SPGrid can be used through the preprocessor magic. The following table shows how
to select different variants of SPGrid:

| GRIDTYPE                  | Refinement                     |
| :------------------------ | :----------------------------- |
| SPGRID                    | Default (Isotropic)            |
| SPGRID_SERIAL             | Default (Isotropic, no MPI)    |
| SPGRID_ISOTROPIC          | Isotropic                      |
| SPGRID_ISOTROPIC_SERIAL   | Isotropic (no MPI)             |
| SPGRID_ANISOTROPIC        | Anisotropic                    |
| SPGRID_ANISOTROPIC_SERIAL | Anisotropic (no MPI)           |
| SPGRID_BISECTION          | Bisection                      |
| SPGRID_BISECTION_SERIAL   | Bisection (no MPI)             |
| SPGRID_COUNT_FLOPS        | use Dune::Fem::Double as ctype |
