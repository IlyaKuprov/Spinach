# experiments/spen/spendosycosy.m

- Signature: `fid=spendosycosy(spin_system,parameters,H,R,K,G,F)`

## Purpose

Simulates the ultrafast 3D DOSY-COSY pulse sequence and returns the free induction decay over both acquisition dimensions and the loop index.

## Physical / mathematical content

- Combines diffusion encoding with COSY coherence transfer. The source applies chirp pulses and gradient intervals, selects the specified coherence orders, and includes the intervening diffusion evolution.
- Uses the Hamiltonian, relaxation, kinetics, gradient, and diffusion/flow superoperators supplied by the imaging context; the measured signal is formed with the supplied detection state.

## Numerical / algorithmic content

- Applies shaped chirp pulses by piecewise-constant propagation and builds propagators for the acquisition periods.
- Stores loop starting states, propagates and detects the signal over both acquired point dimensions, and uses parallel loop execution; GPU arrays are used when enabled.

## Parameters / inputs

- parameters.dims size of the sample in m
- parameters.npts number of spin packets
- parameters.spins nuclei on which the sequence runs
- parameters.deltat timestep for acquisition
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout
- parameters.Ga acquisition gradient in T/m
- parameters.pulsenpoints number of points in the pulse shape
- parameters.smfactor smoothing factor for the pulse
- parameters.Te duration of the pulse
- parameters.Tau duration extra dephasing gradient
- parameters.td diffusion delay, at least
- parameters.Tau+parameters.Te
- parameters.BW bandwidth of the pulse
- parameters.Ge encoding gradient in T/m
- parameters.Gp coherence selection gradient in T/m
- parameters.Tp duration of the coherence selection gradient
- parameters.chirptype can be 'wurst' or 'smoothed'
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow superoperator

## Outputs

- fid -SPENDOSYCOSY free induction decay as a 3D array
- Note: the last five parameters are built automatically by the imaging
- context function.

## Implementation structure

- Forms the Liouvillian and pulse operators, then executes the source-defined diffusion-encoding and COSY preparation with coherence-order selection.
- Builds the acquisition propagators and loop states, and returns an FID array indexed by `npoints1`, `npoints2`, and `nloops`.
