# experiments/spen/spencosy.m

- Signature: `fid=spencosy(spin_system,parameters,H,R,K,G,F)`

## Purpose

Ultrafast COSY pulse sequence with spatial encoding and gradient readout.

## Parameters

- `parameters.dims`: sample size in metres; `parameters.npts`: number of spin packets; `parameters.spins`: nuclei on which the sequence runs.
- `parameters.rho0`: initial state; `parameters.coil`: detection state.
- `parameters.deltat`: acquisition timestep; `parameters.npoints`: acquired points per gradient readout; `parameters.nloops`: loops, each comprising a positive and a negative readout.
- `parameters.Ga`: acquisition gradient in T/m; `parameters.Ge`: encoding gradient in T/m.
- `parameters.pulsenpoints`: pulse-shape points; `parameters.nWURST`: pulse smoothing factor; `parameters.Te`: pulse duration; `parameters.BW`: pulse bandwidth.
- `parameters.Gp`: coherence-selection gradient in T/m; `parameters.Tp`: its duration.
- `parameters.D`: diffusion constant in `m^2/s` (listed in the source header; not accessed directly in this function).
- `H`, `R`, `K`, `G`, `F`: Fokker–Planck Hamiltonian, relaxation, kinetics, gradient, and diffusion/flow operators, respectively. These last five inputs are built automatically by the imaging context function.

## Sequence and output

The function forms `L=H+F+1i*R+1i*K` and applies an initial `pi/2` pulse. It then applies two WURST chirp pulses with opposite encoding gradients, followed by a second `pi/2` pulse bracketed by coherence-selection gradients. After prephasing, it propagates alternating-gradient readout loops and detects `coil'*rho` at each point. The output `fid` is the UFCOSY free induction decay, an array of size `parameters.npoints` by `parameters.nloops`. Loop bodies run in parallel; propagators and states move to a GPU when GPU execution is enabled.

The function requires `sphten-liouv` formalism. It checks that `H`, `R`, `K`, and `F` are equal-sized matrices and that `G` is a cell array.

Authors: jeannicolas.dumez@cnrs.fr; ilya.kuprov@weizmann.ac.il; ludmilla.guduff@cnrs.fr. [Spinach documentation](https://spindynamics.org/wiki/index.php?title=spencosy.m).