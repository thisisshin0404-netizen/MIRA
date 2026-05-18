# MIRA
**M**ultiscale S**I**mulation in **R**adical **A**gony

Multiscale smilation for investigating physical behavior of particles, polymers and biomolecules under extreme conditions like space environment (Not developed yet)

## Build and Run
```bash
cmake -S . -B build
cmake --build build
PYTHONPATH=build python run/test.py
```

## Current MIRA (18 May 2026)
MIRA now supports C++-side initialization of particle positions. The number of particles, box size, and initialization rule are given from Python, while the actual initial placement is generated in C++.

Two initial placement rules are currently available:
- `square`: particles are placed on a square grid fitted to the simulation box
- `random`: particles are placed randomly inside the simulation box

The simulation step still uses periodic boundary conditions, Lennard-Jones 12-6 interactions, and Euler time integration.

## Previous MIRA (12 May 2026)
MIRA was a minimal 2D molecular dynamics toy code written in C++ and exposed to Python with pybind11. It implemented a very basic particle simulation with periodic boundary conditions, Lennard-Jones 12-6 interactions, and Euler time integration.

## Next MIRA
### Atomistic scale
- improve molecular dynamics with a more stable integration algorithm
- generate trajectories
- evaluate the simulation by plotting the Boltzmann distribution

### Coarse-grained / mesoscopic scale
- particle-based coarse-grained modeling is under discussion

### Continuum scale
- PDE-based continuum modeling is under discussion

## Developers

Sarah G. and Gyuhwang S.  
M.Sc. students in Advanced Materials and Processes  
Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)
