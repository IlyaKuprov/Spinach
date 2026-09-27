# examples/spin_chemistry/singlet_yield_3.m

- Signature: `singlet_yield_3()`

## Purpose

Calculate Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore ([DOI](http://dx.doi.org/10.1080/00268979809483134)). The source lists the calculation time as seconds.

## Physical / mathematical content

- Simulates singlet recombination yield for a system of two electrons and one proton (`{'E','E','1H'}`) over a magnetic-field sweep and multiple rate values.
- Sets the electron scalar Zeeman parameters to `2.0023` each and the proton value to `0`. The electron–proton scalar coupling between spins 1 and 3 is `gauss2mhz(20)*1e6`; the spin-3 self-coupling entry is set to `0`.

## Numerical / algorithmic content

- Uses a unit magnet (`sys.magnet=1`) and the `zeeman-hilb` basis with `approximation='none'`.
- Sets `parameters.fields=linspace(0,3*20/1e4,200)` and `parameters.rates=2*pi*[0.005 0.02 0.05 0.1 0.15 0.2 0.3 0.5 2.0]*gauss2mhz(20)*1e6`.
- Sets `parameters.electrons=[1 2]`, `parameters.spins={'E'}`, and `parameters.needs={'zeeman_op'}`. Runs `liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Defines the spin system, basis, couplings, and sequence parameters; then calls `create(sys,inter)` and `basis(spin_system,bas)`.
- Plots the simulated result against `linspace(0,3,200)`, labels the axes `singlet recombination yield` and `$\omega/a$`, and sets axis limits to `[0 3 0.2 1.0]`.
