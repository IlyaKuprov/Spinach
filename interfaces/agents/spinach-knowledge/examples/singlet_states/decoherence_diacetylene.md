# examples/singlet_states/decoherence_diacetylene.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/decoherence_diacetylene.m)

This calculation models diacetylene with two protons and four carbons (4,096-dimensional Liouville space). The spin data are imported from vacuum-DFT coordinates, shifts, couplings and CSAs; the source selects 1H and 13C and passes `[31.8 182.4]` to `g2spinach` as conversion arguments, without assigning units to those values. All dipolar couplings and CSA tensors enter the Redfield relaxation superoperator.

The executable field assignment is 14.1 T, although its immediately preceding comment says 1.0 Tesla; this description follows the assignment. The model sets zero equilibrium, keeps relaxation in the lab frame, uses a 100 ps correlation time and sets both relaxation tolerances to 1e-5. It uses the complete, unapproximated `sphten-liouv` basis.

The function displays 20 small-magnitude relaxation eigenvalues, then constructs the normalised singlet operator for the two centre carbons (spin indices 1 and 2) and evaluates `S'*R*S` as its self-relaxation rate. It also finds two low-magnitude eigenvectors and prints their spherical-tensor composition with `stateinfo`. These are model diagnostics; the source contains no reported numerical result, preparation pulse, gradient, storage-time trace or image reconstruction.
