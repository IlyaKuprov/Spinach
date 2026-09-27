# examples/spin_chemistry/singlet_yield_2.m

- Signature: `singlet_yield_2()`

## Purpose

Simulates a liquid-state magnetic-field effect on a radical pair with six equivalent nuclei using an exponential recombination-kinetics model and full S6 symmetry. The source states a calculation time of seconds.

## Physical / mathematical content

- The system contains two electrons (`E`) and six protons (`1H`). Both electron Zeeman scalar values are `2.002`; the proton values are zero.
- The coupling matrix connects the first electron to each of the six protons. Each listed coupling is `0.295` before the matrix is divided by `2` and converted with `mt2hz`; all other listed entries are zero.
- The kinetics rates are `[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6`. The field array is `1e-3*(0:0.01:5)`.

## Numerical / algorithmic content

- Sets `sys.magnet=1` for the field sweep.
- Uses the `sphten-liouv` formalism with approximation `none`, projections `{0}`, and S6 permutation symmetry for spins 3–8.
- Specifies electrons `[1 2]`, spins `{'E'}`, and the required operator `{'zeeman_op'}`.
- Creates the spin system, applies the basis, and computes `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Defines the system, coupling matrix, symmetry-adapted basis, fields, and kinetics rates.
- Runs the liquid-state simulation and plots `M` against `parameters.fields` as the singlet recombination yield versus magnetic field in Tesla.
