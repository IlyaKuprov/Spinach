# experiments/spen/spendosy.m

- Signature: `fid=spendosy(spin_system,parameters,H,R,K,G,F)`

## Purpose

Simulates the ultrafast DOSY pulse sequence and returns the acquired free induction decay (FID) for each loop.

## Physical / mathematical content

- Forms the evolution generator `L=H+F+1i*R+1i*K`. The sequence uses chirp pulses with the encoding gradient `Ge*G{1}`, selected coherence orders, and positive-gradient intervals of duration `Tau`; the intervening diffusion evolution lasts `td-Tau-Te`.
- During acquisition, the code uses the gradient superoperator with opposite `Ga` polarities for the two readout periods and detects with `parameters.coil`.

## Numerical / algorithmic content

- Applies the chirp waveform with `shaped_pulse_xy` using the `expv-pwc` method, then constructs propagators for the two acquisition-gradient periods.
- Stores the state at each loop start and records `coil'*rho` at each acquired point. Loop bodies run with `parfor`; the state and propagators are moved to the GPU when enabled.

## Parameters / inputs

- parameters.dims size of the sample in m
- parameters.npts number of spin packets
- parameters.spins nuclei on which the sequence runs
- parameters.deltat timestep for acquisition
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout
- parameters.offset offset
- parameters.cond bondary conditions
- parameters.Ga acquisition gradient in T/m
- parameters.pulsenpoints number of points in the pulse shape
- parameters.smfactor smoothing factor for the pulse
- parameters.Te duration of the pulse
- parameters.Tau duration extra dephasing gradient
- parameters.BW bandwidth of the pulse
- parameters.Ge encoding gradient in T/m
- parameters.chirptype can be 'wurst' or 'smoothed'
- parameters.td diffusion delay, at least
- parameters.Tau+parameters.Te
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow superoperator

## Outputs

- fid -free induction decay
- Note: the last five parameters are built automatically by the imaging
- context function.

## Implementation structure

- Checks the Spinach formalism and dimensions, then forms the Liouvillian and the chirp-pulse operators.
- Runs the coherence-selection, encoding, and diffusion-preparation sequence; builds loop propagators and states; and returns an FID array with `npoints` rows and `nloops` columns.
