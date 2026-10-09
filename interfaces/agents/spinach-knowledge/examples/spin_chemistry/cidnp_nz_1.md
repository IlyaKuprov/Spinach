# examples/spin_chemistry/cidnp_nz_1.m

- Signature: `cidnp_nz_1()`
- Source: [examples/spin_chemistry/cidnp_nz_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/cidnp_nz_1.m)

## Purpose

Compares the field dependence of geminate CIDNP from a singlet-born radical pair using Redfield relaxation and a lifetime-shifted Nakajima–Zwanzig kernel. The source comment points to field-cycling CIDNP work at [DOI 10.1039/C3CP44098K](https://doi.org/10.1039/C3CP44098K) as motivation; that pointer is not evidence that the simulated curve is an experimental measurement.

## Model and calculation

- The spin system is two electrons and one proton, `{'E','E','1H'}`. The initial reactant state is the electron singlet `singlet(spin_system,1,2)`; the proton is included in the coupled spin system.
- The field grid is `[0.002 0.005 0.01 0.02 0.035 0.05 0.075]` T (2–75 mT). The source sets the singlet recombination rate to `3e8` Hz and the rotational correlation time to `1e-9` s (1 ns).
- Electron Zeeman factors are `2.0023` and `2.0034`. The proton has an anisotropic hyperfine tensor entered through `mt2hz([0.6 0.6 3.6])`, with zero Euler angles. Recombination uses an explicit first-order singlet-selector loss record at rate `k_rec`, giving the Haberkorn singlet drain; no triplet drain is present.
- For each field, the source constructs the Hamiltonian, relaxation and kinetics superoperators, then compares `redfield` with `naka-zwan`. For the latter it explicitly sets `nz_shift=k_rec/2` (the legacy half-sum Haberkorn scalar approximation) and `nz_onshell=false`. It propagates the singlet in a doubled reactant/product space for `200e-9` s; only the reactant block has Hamiltonian and relaxation dynamics, and the kinetics operator transfers population out of the reactants.

## Observable and plot

The reported quantity is the proton `Lz` expectation in the product block, `real(Nz'*rho_prod)`, labelled product nuclear polarisation. The script prints both theory columns against field and plots them versus field in mT, with Redfield and Nakajima–Zwanzig legend entries. It reports calculation time as minutes. The source comment's interpretation that lifetime broadening changes the product polarisation is explanatory context, not a separately established measurement.
