# examples/spin_chemistry/singlet_yield_5.m

- Signature: `singlet_yield_5()`

## Purpose

Reproduce Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: http://dx.doi.org/10.1080/00268979809483134. The source states a calculation time of seconds.

Source comments list ilya.kuprov@weizmann.ac.il, h.j.hogben@chem.ox.ac.uk, and peter.hore@chem.ox.ac.uk.

## Physical / mathematical content

- Specifies two electrons and two `1H` nuclei, with electron Zeeman scalar values of `2.0023` and nuclear values of `0`.
- Sets couplings `(1,3)` and `(2,4)` to `gauss2mhz(20)*1e6`; explicitly sets `(4,4)` to `0`.
- Computes a singlet recombination yield over a magnetic-field sweep for seven kinetic rates.

## Numerical / algorithmic content

- Sets `sys.magnet=1` and uses a `zeeman-hilb` basis with `bas.approximation='none'`.
- Sets `parameters.fields=linspace(0,0.5*20/1e4,200)` and `parameters.rates=2*pi*[1e-4 1e-3 0.01 0.02 0.05 0.1 0.2]*gauss2mhz(20)*1e6`.
- Sets `parameters.electrons=[1 2]`, `parameters.spins={'E'}`, and `parameters.needs={'zeeman_op'}`. Creates the spin system, applies the basis, and runs `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Specify the spin system, basis, magnetic fields, and kinetics.
- Run the `rydmr_exp` simulation with `liquid` in the lab frame.
- Plot `M` against `linspace(0,0.5,200)`, label the axes `singlet recombination yield` and `$\omega/a$`, and set `axis([0 0.5 0.46 0.58])`.
