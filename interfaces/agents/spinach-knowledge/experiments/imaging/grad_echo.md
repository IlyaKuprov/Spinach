# experiments/imaging/grad_echo.m

- Signature: `fid=grad_echo(spin_system,parameters,H,R,K,G,F)`

## Purpose

Gradient echo pulse sequence. Syntax: fid=grad_echo(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.g_amp -the amplitude of gradient T/m
- parameters.g_step_dur -time step duration
- parameters.g_n_steps -number of time steps

## Outputs

- fid -the time domain echo signal

## Implementation structure

- Gradient echo pulse sequence. Syntax:
- fid=grad_echo(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.g_amp -the amplitude of gradient T/m
- parameters.g_step_dur -time step duration
- parameters.g_n_steps -number of time steps
- fid -the time domain echo signal
- Check consistency
- Assemble the Liouvillian
- Make pulse operators
- Hard 90-degree pulse
