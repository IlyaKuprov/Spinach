# examples/dnp_liq/jdnp/energy_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/energy_levels.m`
- Signature: `energy_levels()`
- Total lines: 52

## Purpose

Energy level diagram transition from the Zeeman limit to the exchange coupling limit in a two-electron system.

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Energy level diagram transition from the Zeeman limit to
- the exchange coupling limit in a two-electron system.
- 600 MHz magnet
- Two electrons
- Exaggerate g-factor difference
- Hilbert space calculation
- Spinach housekeeping
- Relevant operators
- omega_j scan from 0 to omega_e
- Get the energy levels
- Diagonalise the Hamiltonian
- Sort and record energies

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `spin()`, `omega_j()`, `energies()`, `kfigure()`, `kxlabel()`, `kylabel()`.
