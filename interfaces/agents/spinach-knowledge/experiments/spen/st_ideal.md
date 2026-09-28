# experiments/spen/st_ideal.m

- Signature: `inten=st_ideal(spin_system,parameters,H,R,K,G,F)`

## Purpose

Computes the signal for the ideal Stejskal-Tanner diffusion-encoding sequence in the notation of Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 with no gaps between events, and returns the absolute first FID point.

## Physical / mathematical content

- This is the ideal Stejskal-Tanner sequence in the notation of Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 The evolution generator is `L=H+F+1i*R+1i*K`; the two gradient periods use `L+g_amp*G{1}` and each last `delta_sml`. The two intervening delays each last `(delta_big-delta_sml)/2`.
- Applies the source-defined excitation and refocusing pulses and detects the resulting state with the supplied coil.

## Numerical / algorithmic content

- Propagates the state through the pulse sequence, gradient intervals, and delays, then evaluates the detected signal. The function returns a scalar absolute signal rather than a sampled time-domain trace.

## Required inputs

Call from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`. The `parameters` structure must contain:

- `rho0` and `coil`: initial and detection states; `spins`: the working spin.
- `npts`: number of spatial grid points.
- `g_amp`: diffusion-gradient amplitude in T/m.
- `delta_sml` and `delta_big`: the small and big Stejskal–Tanner time intervals, respectively, in seconds.

## Outputs

- inten -the absolute value of the first point in
- the free induction decay; this number is
- proportional to the integral of the real
- part of the correctly phased spectrum

## Implementation structure

- Checks the formalism, input operator dimensions, gradient container, spin selection, and required scalar parameters.
- Constructs the Liouvillian, applies the ideal pulse and diffusion-gradient timing, and returns the absolute detected signal.
