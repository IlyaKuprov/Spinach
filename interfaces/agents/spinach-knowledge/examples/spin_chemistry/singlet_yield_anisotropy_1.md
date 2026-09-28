# examples/spin_chemistry/singlet_yield_anisotropy_1.m

- Signature: `singlet_yield_anisotropy_1()`

## Purpose

Calculate singlet yield anisotropy for a radical pair using an exponential recombination kinetics model. The source states a calculation time of seconds.

## Physical / mathematical content

- The spin system contains two electrons and one `1H` nucleus. The electrons have scalar Zeeman values of `2.0023`; the proton's scalar Zeeman value is `0`.
- An electron–proton coupling is specified between spins 1 and 3 as the diagonal matrix `gauss2mhz([5e6 0 0; 0 4e6 0; 0 0 10e6])`.
- The simulation uses a field of `50e-6`, a recombination rate of `2e6`, and electrons `[1 2]`.

## Numerical / algorithmic content

- Set `sys.magnet=1`, then create the spin system in the `zeeman-hilb` formalism with approximation `none`.
- Use the `leb_2ang_rank_71` orientation grid, `spins={'E'}`, `needs={'zeeman_op'}`, and `sum_up=0`.
- Run `[yield,grid]=powder(spin_system,@rydmr_exp,parameters,'labframe')`. Convert the yield cells to a matrix and subtract `sum(yield.*grid.weights)` from the yield.

## Implementation structure

- Build a hull from `grid.betas` and `grid.gammas`. Map the processed yield to Cartesian coordinates using `x=yield.*sin(grid.betas).*cos(grid.gammas)`, `y=yield.*sin(grid.betas).*sin(grid.gammas)`, and `z=yield.*cos(grid.betas)`.
- Colour the triangulated surface by `sqrt(x.^2+y.^2+z.^2)` and plot it with `trisurf(...,'EdgeAlpha',0.25)`, followed by grid, axis, and box formatting.