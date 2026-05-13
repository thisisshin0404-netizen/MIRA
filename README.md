# MIRA
**M**odeling **I**n **R**adical **A**gony

Multiscale smilation for investigating physical behavior of particles and molecules/polymers under extreme conditions like space environment (Not developed yet)

## Current MIRA (12 May 2026)
MIRA is a minimal 2D molecular dynamics toy code written in C++ and exposed to Python with pybind11. It currently implements a very basic particle simulation with periodic boundary conditions, Lennard-Jones 12-6 interactions, and Euler time integration. The project is built with a small CMake-based structure so the module can be imported directly from Python for testing.

## Next MIRA
### Atomistic scale
- improve molecular dynamics (MD) with a more stable integration algorithm
- generate trajectories
- evaluate the simulation against the Boltzmann distribution

### Coarse-grained / mesoscopic scale
- particle-based coarse-grained modeling is under discussion

### Continuum scale
- PDE-based continuum modeling is under discussion

## Developers

Sarah G. and Gyuhwang S.  
M.Sc. students in Advanced Materials and Processes  
Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)
