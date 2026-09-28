# examples/giant_spin/nuclear_relaxation_1.m

- Signature: `nuclear_relaxation_1()`

## Purpose

Calculate proton relaxation rates and a frequency shift for a rapidly relaxing Dy(III) ion using adiabatic elimination. The example uses a specified ligand field; the source comments give a calculation time of minutes.

## Physical model

- The spin system contains an `E16` Dy(III) electron and a `1H` nucleus at a magnetic field of `14.1`. Their Cartesian coordinates are `[0.00 0.00 0.00]` and `[0.00 5.00 7.00]`, respectively.
- The electron g-tensor is constructed as `V'*diag(D)*V`, with principal values `D=[1.325781 1.322640 1.317917]`; the nuclear shift tensor is zero. Spin–orbit corrections to dipolar couplings are enabled with `sys.enable={'sodd'}`.
- MOLCAS ligand-field coefficients of ranks 2, 4, and 6 are converted with `icm2hz` and `stev2sph`, then rotated with two sets of Euler angles obtained from the supplied direction-cosine matrices. The resulting spherical tensors are assigned to `inter.giant.coeff`.
- Separate electron relaxation times are set to `T1e = T2e = 50 fs`; the corresponding proton rates are zero. Relaxation is kept in the laboratory frame with zero equilibrium.

## Calculation

The calculation uses the `sphten-liouv` formalism without basis approximation. It partitions the 1024 Liouville-space states into four pure-nuclear slow states (`1:4`) and electron-involving fast states (`5:1024`). For each orientation in the loaded `leb_2ang_rank_11.mat` grid, it combines the Hamiltonian with non-interacting relaxation, applies `adelim` to eliminate the fast subspace, and adds the weighted result to a `4×4` nuclear relaxation matrix. The orientation loop uses MATLAB `parfor`.

Normalized proton `Lz` and `L+` states in the slow subspace are used to report `1H R1` and `1H R2` from the real parts of their matrix projections, and `1H DFS` from the imaginary part of the `L+` projection. All three displayed values are labelled `Hz`.
