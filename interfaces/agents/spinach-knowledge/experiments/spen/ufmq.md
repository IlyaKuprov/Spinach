# experiments/spen/ufmq.m

- Signature: `fid=ufmq(spin_system,parameters,H,R,K,G,F)`

## Purpose

Simulates the ultrafast multiple-quantum NMR sequence as a literal implementation of Figure 1A in (http://dx.doi.org/10.1002/cphc.201800667), and returns the acquired free induction decay (FID).

## Physical / mathematical content

- The pulse sequence follows Figure 1A in (http://dx.doi.org/10.1002/cphc.201800667). Forms `L=H+F+1i*R+1i*K` and prepares the state with the source-defined 90-degree, delay, and 180-degree pulse sequence. The final preparation pulse depends on whether the requested multiple-quantum order is even or odd.
- Selects the requested coherence order, applies chirp pulses with opposite encoding-gradient polarities, then selects single-quantum coherence for acquisition.

## Numerical / algorithmic content

- Applies the shaped chirp pulses by piecewise-constant propagation and builds propagators for the alternating acquisition-gradient periods.
- Stores each loop's initial state and records the coil-detected signal point by point. Loop bodies use `parfor`; GPU arrays are used when enabled.

## Required inputs

Call from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`; the spin system must use `sphten-liouv`. The `parameters` structure must contain:

- `rho0` and `coil`: initial and detection states; `spins`: a one-element cell array naming the active nucleus.
- `dims`: sample length in metres; `npts`: number of spatial grid points.
- `npoints`: acquired points per gradient readout; `nloops`: readout loops, each with a positive and a negative readout.
- `Ga` and `Ge`: acquisition and encoding gradient amplitudes in T/m; `deltat`: acquisition time step in seconds.
- `pulsenpoints`: points in the chirp shape; `Te`: chirp duration in seconds; `BW`: chirp bandwidth in Hz; `nWURST`: smoothing parameter; `chirptype`: `wurst` or `smoothed`.
- `delay`: interpulse delay in seconds; `mqorder`: selected multiple-quantum coherence order.

## Outputs

- fid -free induction decay of the ultrafast NMR spectrum.

## Implementation structure

- Checks the formalism and required inputs, then constructs the Liouvillian and chirp waveforms.
- Executes the pulse preparation and multiple-quantum coherence selection, applies the gradient-encoded chirps, and generates an FID array with `npoints` rows and `nloops` columns.
