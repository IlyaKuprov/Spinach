# examples/giant_spin/quartet_magn.m

- Signature: `quartet_magn()`

## Purpose

Calculate and plot the sample magnetisation during a finite-speed magnetic-field sweep for a spin-3/2 particle with zero-field splitting. The source estimates a calculation time of seconds.

## Physical model

- The particle is specified as `sys.isotopes={'E4'}`, with an isotropic Zeeman tensor `diag([2 2 2])`.
- The zero-field-splitting parameters are `D=icm2hz(-0.5)` and `E=0.3*D`. The coupling matrix is `zfs2mat(D,E,0,0,0)`.
- The temperature is `inter.temperature=1.0`.

## Calculation and output

- Set `sys.magnet=1.0` Tesla, as required by the source. Create the spin system with `create(sys,inter)` and set the basis using `bas.approximation='none'` and `bas.formalism='zeeman-hilb'`.
- Scan `parameters.fields=[0 1]` with `parameters.npoints=1000` over `parameters.sweep_time=1e-9` seconds. Set `parameters.orientation=[0 0 0]` and `parameters.nstates=4`.
- Call `[fields,z_magn]=fieldscan_magn(spin_system,parameters)` and plot `z_magn` against `fields`, labelling the axes “Magnetic field, Tesla” and “Sample magnetisation”.
