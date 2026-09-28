# examples/spin_chemistry/singlet_yield_1.m

- Signature: `singlet_yield_1()`

## Purpose

Simulate the liquid-state magnetic-field effect on singlet recombination yield for a radical pair with four nuclei, using an exponential recombination kinetics model. The source estimates a calculation time of seconds.

## Physical / mathematical content

- The system contains two electrons and four `1H` nuclei: `{'E','E','1H','1H','1H','1H'}`. Both electron Zeeman scalar values are `2.002`; the nuclear entries are zero.
- The scalar coupling matrix is constructed as `num2cell(mt2hz([...]/2))`. Its nonzero symmetric entries couple electron 1 to nuclei 3 and 4 with values `0.195`, and electron 2 to nuclei 5 and 6 with values `-1.3` and `0.2`, respectively, before the division by two and `mt2hz` conversion.
- The kinetics rates are `[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6`. The field sweep is `1e-3*(0:0.01:5)`, and the electron indices are `[1 2]`.

## Numerical / algorithmic content

- Set `sys.magnet=1` for the field sweep. Use the `sphten-liouv` basis with `bas.approximation='none'`, `bas.projections={0}`, and `S2` permutation symmetry for spins `[3 4]`.
- Set `parameters.spins={'E'}` and request `parameters.needs={'zeeman_op'}`. Disable ZTE with `sys.disable={'zte'}`.
- Create the spin system, apply the basis, and calculate `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

1. Define the magnet setting, isotopes, Zeeman values, and scalar couplings.
2. Configure the basis, field sweep, kinetics rates, and electron indices.
3. Disable ZTE, then call `create` and `basis`.
4. Run `liquid` with `@rydmr_exp` in the lab frame.
5. Plot `M` against `parameters.fields` as a red line, with axes labelled `magnetic field, Tesla` and `singlet recombination yield`.
