# examples/spin_chemistry/singlet_yield_anisotropy_2.m

- Signature: `singlet_yield_anisotropy_2()`

## Purpose

Calculate singlet yield anisotropy for a radical pair using an exponential recombination kinetics model. The source states a calculation time of seconds.

## Spin system and interactions

- Set the unit magnet to `sys.magnet=1` and the isotopes to `{'E','E','14N','14N','1H'}`.
- Use the Zeeman–Hilbert formalism (`bas.formalism='zeeman-hilb'`) with no basis approximation (`bas.approximation='none'`).
- Define three rotation matrices `R1`, `R2`, and `R3` and interaction-eigenvalue matrices `A1=diag([-1.049,-0.996,13.826])`, `A2=diag([-0.305,-0.222,6.872])`, and `A3=diag([-13.850,-9.372,0.143])`.
- Populate a five-spin coupling matrix at pairs `(1,3)`, `(2,4)`, and `(1,5)` with `1e6*gauss2mhz(Ri*Ai*Ri')` for `i=1,2,3`, respectively.
- Set `inter.zeeman.scalar={2.0023 2.0023 0 0 0}`, then create the spin system and apply the basis.

## Simulation parameters and processing

- Set `parameters.npoints=1`, `parameters.fields=50e-6`, `parameters.rates=50e6`, `parameters.electrons=[1 2]`, `parameters.grid='leb_2ang_rank_35'`, `parameters.spins={'E'}`, `parameters.needs={'zeeman_op'}`, and `parameters.sum_up=0`.
- Compute `[yield,grid]=powder(spin_system,@rydmr_exp,parameters,'labframe')`.
- Convert the yield cells to a matrix and subtract `sum(yield.*grid.weights)` from the yield.
- Obtain a hull from `grid.betas` and `grid.gammas`. Form Cartesian plot coordinates by multiplying the processed yield by the corresponding spherical-direction components, then render a colour-mapped `trisurf` with `EdgeAlpha` set to `0.25`; enable `kgrid` and set the axis and box options.