# examples/singlet_states/decoherence_naphthalenetetrone.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/decoherence_naphthalenetetrone.m)

The naphthalenetetrone example builds a four-proton model (256-dimensional Liouville space) from vacuum-DFT coordinates, chemical shifts, J-couplings and CSA tensors. The relaxation superoperator includes every dipolar coupling and CSA tensor. The import selects 1H and passes `31.8` as a conversion argument; the source does not state a unit for that value.

At 1.0 T, the calculation uses Redfield relaxation with zero equilibrium, lab-frame retention and a 100 ps correlation time. It constructs the full `sphten-liouv` basis without approximation and sets the relaxation-integration and zero tolerances to 1e-5. The function then prints 20 small-magnitude relaxation eigenvalues as rates in Hz. It does not explicitly prepare or propagate a singlet state, nor specify an RF pulse, gradient, storage sequence or imaging step; no calculated eigenvalues or lifetime are asserted here.
