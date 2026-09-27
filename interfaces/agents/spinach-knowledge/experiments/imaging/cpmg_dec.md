# experiments/imaging/cpmg_dec.m

- Signature: `mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)`

## Purpose

Runs a Carr–Purcell–Meiboom–Gill (CPMG) sequence on an MRI phantom from the `imaging()` context and returns an image of the specified spin state.

## Physical / mathematical content

The sequence applies an initial `pi/2` rotation about `Ly`, followed by `parameters.npulses` `pi` rotations about `Lx`. Evolution under `B=H+F+1i*R+1i*K` spans `parameters.dec_time`, with half-delays before the first and after the last `pi` pulse.

## Numerical / algorithmic content

The code constructs spatially replicated `Lx` and `Ly` operators and propagates `parameters.rho0` through the pulses and delays using `step()`. Each pulse-to-pulse delay is `parameters.dec_time/parameters.npulses`.

## Outputs

- mri -amplitude of the detection state at each point of the
- sample
- Note: the spin state to be observed should be specified in
- parameters.coil_st, the coil phantom is ignored.

## Implementation structure

After `grumble()` validates the inputs, the function assembles `B`, constructs pulse operators, runs the CPMG echo train, and converts the final state to an image with `fpl2phan(rho,parameters.coil_st{1},parameters.npts)`. `G` is validated but is not otherwise used.
