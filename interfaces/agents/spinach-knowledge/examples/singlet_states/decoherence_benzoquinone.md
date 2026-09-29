# examples/singlet_states/decoherence_benzoquinone.m

- Signature: `decoherence_benzoquinone()`

## Purpose

Build a Redfield relaxation superoperator for the four-proton para-benzoquinone spin system (256-dimensional Liouville space) from vacuum-DFT spin data, then request its 20 smallest-magnitude relaxation eigenvalues.

## Spin system and relaxation model

The source imports coordinates, chemical shifts, J couplings and chemical-shift-anisotropy (CSA) tensors from `../standard_systems/benzoquinone.log` using `gparse` and `g2spinach` with mapping `{{'H','1H'}}`, argument `31.8` and final argument `[]` (the source does not describe these argument roles); the source comment states that the relaxation model accounts for every dipolar coupling and every CSA tensor. The magnetic field is set to 1.0 T. Redfield relaxation uses zero equilibrium, lab-frame retention and a 100 ps correlation time. The full `sphten-liouv` basis is used, with relaxation-integration and relaxation-zero tolerances of `1e-5`.

## Observable and limits

The script constructs `R` and evaluates `eigs(R-speye(size(R)), 20, 'SM')+1`, labelled in the source as the twenty smallest relaxation rates in Hz. It reports rates only; it does not select or identify a singlet eigenmode, prepare or store a singlet state, or compute a lifetime, image or gradient-encoded signal.

The source comment gives a calculation time of seconds.

Source: [examples/singlet_states/decoherence_benzoquinone.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/decoherence_benzoquinone.m)
