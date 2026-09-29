# examples/giant_spin/nuclear_relaxation_1.m

- MATLAB implementation: [examples/giant_spin/nuclear_relaxation_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/nuclear_relaxation_1.m)

- Signature: `nuclear_relaxation_1()`
- Source: `examples/giant_spin/nuclear_relaxation_1.m`

## Model and tensors

The model couples a Dy(III) `E16` giant spin to a proton (`1H`) at the source-specified field value `14.1`; the source does not state a unit for this value. The electron Zeeman tensor is assembled from principal values `[1.325781, 1.322640, 1.317917]` and the supplied direction-cosine matrix `V`. The proton shift tensor is the zero matrix. `sys.enable={'sodd'}` enables the source-described spin-orbit corrections to dipolar couplings.

The source supplies Cartesian coordinates `[0.00 0.00 0.00]` and `[0.00 5.00 7.00]` for the two spins; it does not state coordinate units. MOLCAS ligand-field coefficients are supplied at ranks 2, 4, and 6. The code converts each set with `icm2hz` and `stev2sph`, then applies Wigner rotations from the supplied ligand-frame and molecular-frame direction-cosine matrices before assigning the resulting spherical tensors to `inter.giant.coeff`. The input-coefficient units are not stated in the source.

Separate electron relaxation times are set to T1e = T2e = 50 fs (rates `1/50e-15`); the proton relaxation rates are zero. Relaxation is kept in the lab frame with zero equilibrium. The basis uses `sphten-liouv` and `approximation='none'`.

## Orientation average and reported quantities

The 1024 Liouville-space states are partitioned into four pure-nuclear slow states (indices 1:4) and electron-involving fast states (5:1024). For each weighted orientation in `leb_2ang_rank_11.mat`, the example forms `L = I + orientation(Q,angles) + 1i*R_ni`, where `R_ni` is the non-interacting relaxation superoperator, and uses `adelim(spin_system,L,fast_idx,slow_idx)` to obtain the effective slow-subspace relaxation matrix. MATLAB `parfor` performs the orientation loop and the weighted matrices are accumulated.

Normalised proton `Lz` and `L+` states restricted to the slow subspace give the displayed `1H R1` and `1H R2` from the real parts of the negative matrix projections, and `1H DFS` from the imaginary part of the negative `L+` projection. All three printed values are labelled Hz. The source estimates a calculation time of minutes; it contains no fixed numerical output values or plot.
