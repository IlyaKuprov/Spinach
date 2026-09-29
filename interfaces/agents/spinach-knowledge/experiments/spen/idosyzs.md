# experiments/spen/idosyzs.m

- MATLAB implementation: [experiments/spen/idosyzs.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/idosyzs.m)

- Signature: `inten=idosyzs(spin_system,parameters,H,R,K,G,F)`
- Source documentation: <https://spindynamics.org/wiki/index.php?title=idosyzs.m>
- Shipped example: `examples/nmr_spen/idosyzs_test_1.m`

## Purpose and context

This is a simplified model of the Zangger–Sterk (ZS) iDOSY pulse sequence. It is called through `imaging()`, which prepares and supplies the sequence operators `H`, `R`, `K`, `G`, and `F`, together with the spatially encoded initial state and detection coil in `parameters`.

## Parameters

- `rho0`: initial state; `coil`: detection state; `spins`: nuclei on which the sequence operates (the code uses `spins{1}`).
- `g_amp`: diffusion-encoding gradient amplitude (T/m); `sel_g_amp`: gradient amplitude during the selective pulse (T/m).
- `rf_phi`: inversion-pulse flip angle in radians. Although the source parameter comment calls this a “phase” and labels it “rad/s”, the implementation divides it by a duration to obtain an RF amplitude in rad/s. The shipped example sets `rf_phi=pi` for this inversion pulse. It is therefore an angle input, not a phase rate.
- `rf_dur`: duration of the selective 180-degree shaped pulse (s).
- `delta_big`: diffusion delay (s); `delta_sml`: diffusion-encoding gradient pulse width (s).
- `filename`: character-string pulse-shape name; `pulse_npoints`: number of points used to sample the soft pulse.
- `diff`: diffusion coefficient (m^2/s); `dims`: sample dimension; `npts`: number of spin packets. These fields are required by the function's input checks and used in the surrounding imaging setup.

## Sequence and RF calibration

The sequence forms `L=H+F+1i*R+1i*K`, and constructs `Lx` and `Ly` from the raising operator for `spins{1}`, expanded over the spin-packet dimensions. It sets the free diffusion interval to `delta=delta_big-delta_sml-rf_dur`, reads the waveform with `read_wave(filename,pulse_npoints)`, and divides `rf_dur` into equal time steps.

The RF calibration is `gamma_B1=rf_phi/rf_dur` (rad/s), followed by `gamma_B1_max=gamma_B1/scaling_factor`. The waveform amplitudes are scaled by this value, then its amplitudes and phases are converted to `Cx,Cy`. The reported field is `gamma_B1_max/(2*pi)` in Hz. Thus the duration converts the flip-angle input into a rate; the source comment's “phase ... (rad/s)” wording does not describe the value used by the executable code.

The propagation order is:

1. Apply a `pi/2` pulse using `Ly` to `rho0`, then select coherence `-1` on `spins{1}`.
2. Evolve for `delta_sml` under the first **positive** diffusion gradient, `L+g_amp*G{1}`, then evolve under `L` for `delta`.
3. Apply the shaped 180-degree pulse using `shaped_pulse_xy`, with `{Lx,Ly,G{1}}`, RF components `{Cx,Cy}`, and a positive selective-gradient amplitude `+sel_g_amp` on `G{1}` at each pulse point. The propagator uses `expv-pwc`.
4. Select coherence `+1` on `spins{1}`, apply the second **positive** `L+g_amp*G{1}` diffusion gradient for `delta_sml`, and evolve under `L` for `delta` to refocus the signal.

## Output and checks

The returned signal is `inten=abs(coil'*rho)`, the magnitude of the detected first FID point. In the imaging context this is proportional to the real-spectrum integral for a correctly phased spectrum. The routine reports the waveform scaling factor and the maximum RF field in Hz.

Input checks require `sphten-liouv` formalism; numeric matrix inputs `H`, `R`, `K`, and `F` of equal size; and cell-valued `G`. The fields `rho0`, `coil`, `dims`, `npts`, `spins`, `filename`, `pulse_npoints`, `diff`, `g_amp`, `sel_g_amp`, `rf_dur`, `rf_phi`, `delta_big`, and `delta_sml` must be present. The scalar fields `dims`, `npts`, `pulse_npoints`, `diff`, `g_amp`, `sel_g_amp`, `rf_dur`, `rf_phi`, `delta_big`, and `delta_sml` are each required to contain one element; `spins` must also contain one element; `filename` must be a character string. It also checks `delta_big>delta_sml+rf_dur`.

## Shipped example

`examples/nmr_spen/idosyzs_test_1.m` uses `rf_phi=pi`, `rf_dur=0.045 s`, and `gaussian_1000.pk` sampled at 100 points. It creates a 4000-point phantom and evaluates 20 diffusion-gradient amplitudes from 0.01 to 0.40 T/m. This concrete setup is consistent with `rf_phi` being an inversion flip angle in radians, while the implementation's division by `rf_dur` supplies the RF rate.
