# examples/dnp_sol/cavity_dnp/cavity_dnp_minimal.m

Demonstrates one ideal cavity-assisted DNP cycle for a single Bi:Si donor using three active energy eigenstates and an explicit truncated bosonic cavity. The full donor Hamiltonian identifies the active states; seven inert ground populations are included analytically in normalisation.

Call `result=cavity_dnp_minimal(nlevels)` with an integer Fock dimension of at least two. The returned structure contains time in seconds, active populations, photon occupation, the target-state population enhancement, the change in physical nuclear magnetisation, and density-matrix consistency diagnostics. Frequencies used internally are angular frequencies. No files or figures are produced.

The pulses are ideal selective pi rotations with donor-cavity coupling switched off. The waiting stage retains Jaynes-Cummings coupling and photon loss. The illustrative coupling is not fitted to an ensemble; collective effects and repeated pumping are outside the model. See the adjacent example README for parameters and the reference to Lorriaux et al., arXiv:2610.19180v1.
