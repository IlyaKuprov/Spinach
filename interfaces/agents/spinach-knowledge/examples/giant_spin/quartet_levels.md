# examples/giant_spin/quartet_levels.m

- Signature: `quartet_levels()`

## Purpose

Calculate the four energy levels of a spin-3/2 particle with zero-field splitting as the magnetic field is scanned from 0 to 1 T. The calculation takes seconds.

## Physical / mathematical content

The particle is specified as `E4` with an isotropic Zeeman tensor, `diag([2 2 2])`. The zero-field-splitting parameters are `D=icm2hz(-0.5)` and `E=0.3*D`; `zfs2mat(D,E,0,0,0)` constructs the coupling matrix with all three orientation angles set to zero. `sys.magnet` is set to `1.0` T, as required by the example.

## Numerical / algorithmic content

The calculation uses the `zeeman-hilb` formalism with `bas.approximation='none'`. `fieldscan_enlev` evaluates four energy levels at `parameters.npoints=100` points over `parameters.fields=[0 1]`, with `parameters.orientation=[0 0 0]`.

## Implementation structure

The function `quartet_levels()` defines the particle, Zeeman tensor and zero-field-splitting matrix; creates the Spinach spin system with `create(sys,inter)`; applies the basis with `basis(spin_system,bas)`; and calls `fieldscan_enlev(spin_system,parameters)`.
