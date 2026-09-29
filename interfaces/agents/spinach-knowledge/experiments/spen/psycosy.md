# experiments/spen/psycosy.m

Canonical source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/psycosy.m
Source paper: [Kenwright, spatially encoded COSY, Figure 4](https://doi.org/10.1002/mrc.4727)
Spinach Wiki: https://spindynamics.org/wiki/index.php?title=psycosy.m

`fid=psycosy(spin_system,parameters,H,R,K,G,F)` implements Alan Kenwright's spatially encoded COSY sequence and returns a two-dimensional free-induction decay. The source examples for `spins` include `{'1H'}` and `{'13C'}`. The imaging context supplies `H`, `R`, `K`, `G`, and `F`; the sequence uses `L=H+F+1i*R+1i*K` and spatially extends its `Lx` and `Ly` pulse operators.

## Model and sampling

The model starts from `rho0` with a hard 90-degree pulse about `Lx`. F1 evolution is split into two halves: the first is sampled in trajectory mode, followed by selection of `+1` coherence, a hard 180-degree pulse, and selection of `-1` coherence. The spatial-encoding element is a saltire chirp pair under `L+gamp*G{1}`. Its first chirp uses `{Cx,+Cy}`; after evolution under that same gradient-containing generator for `sal_del-2*sal_dur`, the sequence selects zero coherence and applies the second chirp with `{Cx,-Cy}`. The source generates the saltire waveform from `sal_npt`, `sal_dur`, `sal_swp`, and `sal_smf`, normalises it by `max(Cx)`, and sets the RF scale using `sal_ang`.

The second F1 half uses refocusing mode with timestep `0.5/sweep` and `npoints(1)-1` steps. After selecting `+1` coherence, the model evolves for `tmix`, selects `+1` again, and applies a final hard 90-degree pulse about `Lx`. F2 is acquired as an observable evolution under `L` with `coil`, timestep `1/sweep`, and `npoints(2)-1` steps. The returned FID has F1 and F2 point counts `npoints(1)` and `npoints(2)`; no spectrum or experimental measurement is produced by this function.

The timing inputs are in seconds where applicable: `tmix`, `sal_dur`, and `sal_del`. The interval between the two chirp waveforms is explicitly `sal_del-2*sal_dur`, with the input check `sal_del>2*sal_dur`. `sweep` is in Hz, `sal_swp` is in Hz, `sal_ang` is in degrees, and `gamp` is in T/m. The source also assigns `parameters.delta=1/(4*parameters.sweep)`; its subsequent propagation calls use the explicit F1 and F2 timings described above rather than reading that field.

## Required inputs

The function checks for `rho0`, `coil`, one-element `spins`, `sweep`, two-element `npoints` with integer lengths greater than one, `tmix`, `gamp`, `sal_ang`, `sal_dur`, `sal_del`, `sal_swp`, `sal_npt`, and `sal_smf`. It requires positive `sweep`, `sal_dur`, and `sal_swp`; non-negative `tmix`; `sal_del>2*sal_dur`; positive-integer `sal_npt`; and `sal_smf` from 0 through 50. `gamp` and `sal_ang` must be finite real scalars. The body also reads `parameters.npts` from the imaging-grid context to replicate pulse operators. It requires `sphten-liouv` formalism, same-sized matrix operators `H`, `R`, `K`, and `F`, and a three-element cell array `G` of numeric gradient operators; this sequence uses `G{1}`.
