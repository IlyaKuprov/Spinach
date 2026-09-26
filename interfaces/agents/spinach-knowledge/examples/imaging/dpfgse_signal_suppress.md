# examples/imaging/dpfgse_signal_suppress.m

- Signature: `dpfgse_signal_suppress()`

## Purpose

DPFGSE water suppression example for a solution of GABA in water. Gradients and soft pulses are done explicitly. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- DPFGSE water suppression example for a solution of GABA
- in water. Gradients and soft pulses are done explicitly.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnet (Tesla)
- Chemical shifts (ppm)
- J-couplings (Hz)
- Basis set
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Spinach housekeeping
