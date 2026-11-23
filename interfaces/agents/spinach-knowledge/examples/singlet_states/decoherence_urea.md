# examples/singlet_states/decoherence_urea.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/decoherence_urea.m)

This urea example is framed in the source as a demonstration that its nitrogen singlet is not long-lived. It imports hydrogen and nitrogen from a vacuum-DFT calculation, then uses `isoswap` to replace the imported `14N` with `15N` and rescale its couplings. The Redfield relaxation superoperator includes every dipolar coupling and CSA tensor. The source does not state a spin count or Liouville-space dimension. Its conversion call passes `[30.0 166.0]` for hydrogen and nitrogen; no units are specified there.

The model uses a 1.0 T field, zero equilibrium, lab-frame relaxation, a 100 ps correlation time and 1e-5 integration and zero tolerances, with the complete unapproximated `sphten-liouv` basis. It constructs the 15N longitudinal operator `Lz` and the singlet operator on spin indices 1 and 4, then prints `norm(R*Lz)/norm(Lz)` and `norm(R*S)/norm(S)`. Those expressions describe the diagnostics performed; no output values or lifetime are supplied here. The function contains no RF-preparation pulse, gradient, time-domain storage sequence or imaging reconstruction.
