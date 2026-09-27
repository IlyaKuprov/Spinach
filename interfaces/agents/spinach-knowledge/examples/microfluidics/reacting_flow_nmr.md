# examples/microfluidics/reacting_flow_nmr.m

- Signature: `reacting_flow_nmr()`

## Purpose

This example couples flow and diffusion through a microfluidic chip to two second-order Diels–Alder reaction pathways, then simulates proton NMR detection in the narrow region assumed to contain the coil. The calculation takes days on a CPU and is much faster on a GPU.

## Physical / mathematical content

COMSOL mesh and velocity data define transport across the chip. Cell-wise concentrations of the reactants and products evolve under advection, diffusion, and competing exo/endo reaction channels; the NMR signal is detected only over the selected coil footprint.

## Numerical / algorithmic content

A shared flow-diffusion generator is combined at each time step with local reaction generators to advance the concentration field. The resulting cell concentrations are interpolated in time and used to initialise spin states at successive reaction times. During each NMR acquisition, the spin system evolves with transport and concentration-dependent reaction terms, using a two-point Lie quadrature step; the FIDs are then apodised, zero-filled, Fourier transformed, and displayed as a waterfall spectrum.

## Implementation structure

The function imports the Diels–Alder kinetics and COMSOL hydrodynamics, builds the mesh and spin system, and advances the concentration dynamics before running the coupled chemistry, hydrodynamics, and spin-dynamics stage. A rectangular phantom selects the detection region, and spectra are calculated at sampled reaction times.
