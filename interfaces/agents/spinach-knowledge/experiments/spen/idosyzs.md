# experiments/spen/idosyzs.m

- Signature: `inten=idosyzs(spin_system,parameters,H,R,K,G,F)`
- Source: <https://spindynamics.org/wiki/index.php?title=idosyzs.m>

A simplified ZS iDOSY sequence, called from `imaging()`, which supplies `H`, `R`, `K`, `G`, and `F`.

## Parameters

- `rho0`: initial state; `coil`: detection state; `spins`: nuclei on which the sequence runs.
- `g_amp`: diffusion-encoding gradient amplitude (T/m); `sel_g_amp`: gradient amplitude during the selective pulse (T/m).
- `rf_phi`: phase of the inversion pulse (rad/s), as labeled in the source comment; `rf_dur`: selective 180-degree pulse duration (s).
- `delta_big`: diffusion delay (s); `delta_sml`: diffusion-encoding gradient pulse width (s).
- `filename`: pulse-shape character string; `pulse_npoints`: number of soft-pulse points; `diff`: diffusion constant (m^2/s).
- `dims`: sample dimension; `npts`: number of spin packets.

## Sequence

The code forms `L=H+F+1i*R+1i*K` and constructs `Lx` and `Ly` from the raising operator for `spins{1}`. It sets `delta=delta_big-delta_sml-rf_dur`, reads the waveform with `read_wave(filename,pulse_npoints)`, and divides `rf_dur` into `pulse_npoints` equal time steps. For RF calibration it computes `gamma_B1=rf_phi/rf_dur` (rad/s), `gamma_B1_max=gamma_B1/scaling_factor`, scales the waveform amplitudes, and converts amplitudes and phases to `Cx,Cy`.

1. Apply a `pi/2` pulse using `Ly` to `rho0`, then select coherence `-1` on `spins{1}`.
2. Evolve for `delta_sml` under the **positive** diffusion-encoding gradient `L+g_amp*G{1}`, then evolve under `L` for `delta`.
3. Apply the shaped 180-degree pulse with `shaped_pulse_xy`, using `{Lx,Ly,G{1}}`, `{Cx,Cy,+gradient_amplitudes}`, and `expv-pwc`; every selective-pulse gradient amplitude is `+sel_g_amp` on `G{1}`.
4. Select coherence `+1` on `spins{1}`, evolve for `delta_sml` under the second **positive** `L+g_amp*G{1}` gradient, then evolve under `L` for `delta` to refocus the signal.

## Output and checks

`inten=abs(coil'*rho)` is the absolute first FID point, proportional to the integral of the real part of a correctly phased spectrum. The function reports the waveform scaling factor and `gamma_B1_max/(2*pi)` in Hz.

Validation requires `sphten-liouv` formalism; numeric, matrix-valued `H`, `R`, `K`, and `F` of equal size; and a cell-valued `G`. It requires the listed parameter fields, with exactly one element in each of `dims`, `npts`, `spins`, `pulse_npoints`, `diff`, `g_amp`, `sel_g_amp`, `rf_dur`, `rf_phi`, `delta_big`, and `delta_sml`; `filename` must be a character string. It also requires `delta_big>delta_sml+rf_dur`.