# MIRA
**M**ultiscale S**I**mulation in **R**adical **A**gony

Multiscale smilation for investigating physical behavior of polymers and biomolecules for building parameters and modeling biological phenomena (Not developed yet)

## Build and Run
```bash
cmake -S . -B build
cmake --build build
PYTHONPATH=build python run/test.py
```

## Current MIRA (28 May 2026)
MIRA now supports random C++-side initialization of particle positions. The number of particles, box size, and random seed are given from Python, while the actual initial placement is generated in C++.

The simulation can now be run for a user-defined number of time steps. The number of steps is controlled from Python, while each simulation step is computed in C++.

Two time integration rules are currently available:
- `euler`: semi-implicit Euler time integration
- `verlet`: velocity Verlet time integration

The simulation step uses periodic boundary conditions and Lennard-Jones 12-6 interactions.

## Previous MIRA 
(18 May 2026) MIRA supported C++-side initialization of particle positions. The number of particles, box size, and initialization rule were given from Python, while the actual initial placement was generated in C++. The simulation step used periodic boundary conditions, Lennard-Jones 12-6 interactions, and Euler time integration.

(12 May 2026) MIRA was a minimal 2D molecular dynamics toy code written in C++ and exposed to Python with pybind11.

## Next MIRA
### Atomistic scale
- improve molecular dynamics with a more stable integration algorithm
- generate trajectories
- evaluate the simulation by plotting the Boltzmann distribution

### Mesoscopic scale
- particle-based coarse-grained modeling is under discussion
- Discrete Element Method is under disccusion

## Developers

Sarah G. and Gyuhwang S.  
M.Sc. students in Advanced Materials and Processes  
Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)
