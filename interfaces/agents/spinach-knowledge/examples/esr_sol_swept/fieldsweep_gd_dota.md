# examples/esr_sol_swept/fieldsweep_gd_dota.m

- Signature: `fieldsweep_gd_dota()`

## Purpose

Simulate a powder-averaged, W-band field-swept ESR spectrum of a Gd(III) DOTA complex using exact diagonalisation. The source notes a calculation time of seconds.

## Physical / mathematical content

- Models the electron spin with `E8`, a scalar Zeeman parameter of `1.9918`, and a traceless axial self-coupling tensor with principal values `[0.57e9, 0.57e9, -2*0.57e9]/3`.
- Uses a 90 GHz microwave frequency and sweeps magnetic field from 3.05 to 3.4 T. The source sets `sys.magnet=1` (commented “must be 1”). No hyperfine tensor or anisotropic g tensor is specified.

## Numerical / algorithmic content

- Uses the `zeeman-hilb` formalism with `bas.approximation='none'` for exact diagonalisation.
- Performs powder averaging on the `rep_2ang_100pts_sph` orientation grid. Sets a linewidth of `2e-4` T, 4096 field points, `int_tol=10.0`, `tm_tol=0.1`, and `rspt_order=Inf`.
- Sets the high-temperature initial state to `-state(spin_system,'Lz','E8')` and computes the spectrum with `fieldsweep`.

## Implementation structure

The function creates the `E8` spin system, assigns its Zeeman and self-coupling parameters, builds the basis, runs `fieldsweep`, and plots intensity against magnetic field in tesla.