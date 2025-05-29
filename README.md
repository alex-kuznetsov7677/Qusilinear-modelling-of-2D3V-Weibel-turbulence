# Quasilinear modelling of 2D3V Weibel turbulence
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Publication](https://img.shields.io/badge/Publication-ResearchGate-b31b1b.svg)](https://www.researchgate.net/publication/377742346_Quasilinear_Simulation_of_the_Development_of_Weibel_Turbulence_in_Anisotropic_Collisionless_Plasma)

## 📌 Key Features
- Solves the Vlasov-Maxwell equations in 3D velocity-space and 2D wavevector-space grids.
- Simulates the quasilinear interaction between modes of Weibel turbulence in a weakly collisional, magnetized, anisotropic plasma
- Tracks spectral dynamics and field energy transfer
- Implements Stormer-Verlet (Leapfrog) numerical scheme
- MPI-parallelized for efficient computation of mode dynamics across distributed systems (e.g., KIAM K100 supercomputer).

## Compilation and Execution

This code was developed on the [KIAM K100 supercomputer (Moscow, Russia)](https://www.kiam.ru/MVS/).

  To compile the program, run:
  ```bash
  make
  ```
  After compilation, run the program with:
  ```bash
  mpirun -np <number_of_processes> -maxtime <time_of_execution_in_minutes> build/exe9_c
  ```
## Configuration

Simulation parameters are set in config.txt. Key parameters include:

- start_level_modes  &mdash;  initial perturbation level for each Weibel mode
  
- Ngarmonik_r, Ngarmonik_phi  &mdash;  number of filamentation modes in r and phi in polar coordinate
  
- Beta_perp &mdash;  thermal velocities of background and beam plasma
- A  &mdash;  relative beam density
- MAGNIT  &mdash;  external uniform magnetic field
- type_of_collisions  &mdash; 0 - without collisions; 1 - anomolous collisions proportional to `b²/(8π)`
  
- setkaBB  &mdash;  velocity space grid size
- ratio_Vstepx_multiple_setkaBB_to_Beta_perp  &mdash;  X-step scaling - Determines velocity grid spacing in x-direction: Δv<sub>x</sub> = (β<sub>bkg</sub> × 8)/setkaBB

- Tmax  maximum simulation time (in plasma frequency)
- dt  time step (in plasma frequency)

- kmin_r  &mdash;  Minimum k<sub>r</sub> - Lower bound of wavevector range in r-direction
- kstep_r  &mdash;  k<sub>r</sub> step - Spacing between wavevectors in r-direction

## License

Apache License 2.0 - See [LICENSE](LICENSE) for details.

## ❓ Support

For questions and bug reports:

- [Open an issue](https://github.com/alex-kuznetsov7677/Quasilinear-Weibel-and-Lengmuir-Turbulence)
- Email:
  - [kuznetsov.alexey@ipfran.ru](mailto:kuznetsov.alexey@ipfran.ru)
  - [a.kuznetsov7677@gmail.com](mailto:a.kuznetsov7677@gmail.com)
