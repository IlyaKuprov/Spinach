# experiments/spen/spendosy.m

- Signature: `fid=spendosy(spin_system,parameters,H,R,K,G,F)`

## Purpose

Ultrafast DOSY pulse sequence. Syntax: fid=spendosy(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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

- Ultrafast DOSY pulse sequence. Syntax:
- fid=spendosy(spin_system,parameters,H,R,K,G,F)
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
