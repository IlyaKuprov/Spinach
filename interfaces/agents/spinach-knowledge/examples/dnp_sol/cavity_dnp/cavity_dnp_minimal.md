# examples/dnp_sol/cavity_dnp/cavity_dnp_minimal.m

Illustrates one ideal cavity-assisted DNP cycle for Bi:Si using the native Spinach particle specification `E`, `209Bi`, and `Cn`. The full donor-cavity Hilbert space has dimension `20*nlevels`; all donor spectator levels and retained cavity occupations remain represented.

Call `result=cavity_dnp_minimal(nlevels)` with an integer cavity dimension of at least three. The result contains time in seconds, active-state populations, photon occupation, target-state population enhancement, physical nuclear magnetisation change, diagnostics, and the native system and generators. No files or figures are produced.

`create`, `basis`, and `hamiltonian` identify the physical donor energy eigenstates. Native `operator` and `state` calls provide composite observables and vacuum state. Spectral projectors specify ideal selective microwave/RF transitions; `hilb2liouv`, native mode `relaxation`, and `evolution` provide commutators, photon loss, and propagation. No external reduced density solver or custom tensor dissipator is used.

The cavity remains explicit during emission, but donor-cavity coupling is switched off during ideal pulses. Coupling and pulse strengths are illustrative rather than fitted; collective effects, disorder, repeated pumping, spin diffusion, and non-cavity relaxation are excluded. See the example README and Lorriaux et al., arXiv:2610.19180v1.
