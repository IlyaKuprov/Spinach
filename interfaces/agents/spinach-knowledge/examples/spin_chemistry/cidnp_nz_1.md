# examples/spin_chemistry/cidnp_nz_1.m

- Signature: `cidnp_nz_1()`

## Purpose

Field dependence of geminate CIDNP from a radical pair in a viscous solvent, computed with Redfield theory and the lifetime-shifted Nakajima-Zwanzig kernel. At a rotational correlation time of 1 ns and a singlet recombination rate of `3e8 Hz`, anisotropic hyperfine relaxation proceeds at a fair fraction of the recombination rate, and the pair drains before the bath decorrelates. Lifetime broadening of the spectral densities changes the nuclear polarisation left in the diamagnetic product. The doubled-space setup follows `cidnp_geminate.m`; the low-field regime is motivated by the field-cycling CIDNP work of the Yurkovskaya and Ivanov school (https://doi.org/10.1039/c3cp44098k). Calculation time: minutes

## Physical / mathematical content


## Physical / mathematical content

- The field grid is `[0.002 0.005 0.01 0.02 0.035 0.05 0.075] T`; the singlet recombination rate is `3e8 Hz`, and `tau_c=1e-9 s`. The model has two electrons and one proton, with electron g factors `2.0023` and `2.0034`, and an anisotropic proton hyperfine tensor from `mt2hz([0.6 0.6 3.6])`.
- Both Redfield and Nakajima-Zwanzig calculations retain lab-frame relaxation terms and DFS terms. For the latter, chemical lifetime shifting is enabled off shell.

## Numerical / algorithmic content

- At each field, the script constructs a doubled reactant/product problem, starts from an electron singlet, and propagates the coupled Hamiltonian, relaxation, and reaction dynamics for `200 ns`. It records product-block proton `Lz` polarisation and plots both theory curves against field.

## Implementation structure

- Field dependence of geminate CIDNP from a radical pair in a viscous
- solvent, computed with Redfield theory and with the lifetime-shifted
- Nakajima-Zwanzig kernel. At a rotational correlation time of 1 ns and
- a singlet recombination rate of 3e8 Hz, the anisotropic hyperfine
- relaxation proceeds at a fair fraction of the recombination rate, and
- the pair drains before the bath decorrelates: the spectral densities
- seen by the surviving pair are lifetime-broadened, which changes the
- nuclear polarisation left in the diamagnetic product. The doubled-space
- bookkeeping follows the cidnp_geminate.m example; the low-field regime
- is motivated by the field-cycling CIDNP work of the Yurkovskaya and
- Ivanov school:
- Calculation time: minutes
