# experiments/spen/spencosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/spencosy.m`
- Signature: `fid=spencosy(spin_system,parameters,H,R,K,G,F)`
- Total lines: 245

## Purpose

Ultrafast COSY pulse sequence. Syntax: fid=spencosy(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
- parameters.nWURST smoothing factor for the pulse
- parameters.Te duration of the pulse
- parameters.BW bandwidth of the pulse
- parameters.Ge encoding gradient in T/m
- parameters.Gp coherence selection gradient in T/m
- parameters.Tp duration of the coherence selection gradient
- parameters.D diffusion constant, m^2/s
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow superoperator

## Outputs

- fid UFCOSY free induction decay
- Note: the last five parameters are built automatically by the imaging
- context function.

## Implementation structure

- Ultrafast COSY pulse sequence. Syntax:
- fid=spencosy(spin_system,parameters,H,R,K,G,F)
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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `step()`, `report()`, `shaped_pulse_xy()`, `propagator()`, `clean_up()`, `clear()`, `ismember()`, `gpuArray()`, `local_fid()`, `gather()`, `fid()`, `ismatrix()`.
