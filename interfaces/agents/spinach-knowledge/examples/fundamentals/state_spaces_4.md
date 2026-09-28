# examples/fundamentals/state_spaces_4.m

- Signature: `state_spaces_4()`

## Purpose

Calculate and analyse a magic-angle-spinning (MAS) trajectory for an isotopically labelled glycine powder, starting from proton L+ magnetisation. The source estimates hours of calculation time and notes that a GPU can make it faster.

## Method

The script reads `../standard_systems/glycine.log` with `gparse`/`g2spinach`, specifying 1H, 13C, and 15N and the shift values `[31.5 182.1 264.5]`. It uses a 14.1 T field, the `sphten-liouv` formalism, no basis approximation, longitudinal `{15N,13C}`, and `sys.tols.krylov_tol=1000`. MAS parameters are rate 2000 Hz, axis `[1 1 1]`, `max_rank=17`, sweep `1e5`, 64 points, offset 15000, 13C detection, and grid `rep_2ang_200pts_sph`. The initial state and coil are both proton L+.

`singlerot` generates the trajectory; `fpl2rho` averages over rotor phase before `trajan` plots correlation order. The `sys.enable={'gpu'}` line is commented out, so this script does not enable GPU execution as written.
