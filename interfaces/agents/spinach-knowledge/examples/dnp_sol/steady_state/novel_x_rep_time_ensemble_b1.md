# examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1.m

- Signature: `novel_x_rep_time_ensemble_b1()`
- Source: [`examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1.m)

## Purpose

Compares steady-state NOVEL proton polarisation across a repetition-time scan, with and without a flipback pulse, while averaging over a six-node microwave B1 distribution. The source estimates hours of calculation time. This is a parameterised example script, not a function that accepts scan settings or returns the computed arrays.

## Spin system and relaxation

The example labels the magnet as X-band and sets `sys.magnet=0.34`. It contains an electron and one proton. The electron Zeeman eigenvalues are [2.00319, 2.00319, 2.00258]; the proton shift is [0,0,5] (the source calls this a ppm guess). Euler angles are set to `(pi/180)*{[0,10,0],[0,0,10]}`. Temperature is 80. Coordinates are [0,0,0] for the electron and [0,0,3.500] for the proton; the source derives `r_en` from the latter coordinate but does not annotate coordinate units.

The basis uses `sphten-liouv` without an approximation, and the propagator chop tolerance is `1e-12`. The source disables `hygiene`. Relaxation uses `t1_t2`, diagonal retention, and DiBari equilibrium. Proton R1 is a function handle to `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,26,r_en,bet)`; the electron R1 is 1e3, while the two R2 values are 200e3 and 50e3. The relaxation helper and steady-state callback are not defined in this file: `r1n_dnp` is defined in [`etc/textbook/r1n_dnp.m`](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r1n_dnp.m), and `noveldnp_steady` in [`experiments/hyperpol/noveldnp_steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/noveldnp_steady.m). Both implementations must be on the MATLAB path.

## B1 and repetition-time calculation

The B1 quadrature nodes and weights come from `gaussleg(14e6,16e6,5)`, with the interval explicitly labelled Hz. Repetition times are 30 logarithmically spaced values from 1e-4 to 1e-2; the plot converts them to milliseconds, giving 0.1–10 ms. Each B1 node sets `irr_powers` and a 90-degree pulse duration of `1/(4*irr_powers)`. For every repetition time, the code performs two ESR-context powder calculations on `rep_2ang_800pts_sph`: one without flipback and one with flipback. Both use a 500 ns contact pulse; the shot spacing is set to repetition time minus the relevant pulse/contact duration (one pulse without flipback, two with flipback). The parameters also set `flippulse=1` (commented as NOVEL; 0 denotes solid effect), `addshift=-3.3e6`, and `el_offs=0e6`.

## Result and limits

The two result arrays are integrated over B1 using the quadrature weights normalised by their sum. The figure plots the real proton longitudinal expectation value against repetition time for the no-flipback and flipback cases, and is saved as `novel_x_rep_time_ensemble_b1.fig`. The source specifies the computation and plot but contains no numerical polarisation results.
