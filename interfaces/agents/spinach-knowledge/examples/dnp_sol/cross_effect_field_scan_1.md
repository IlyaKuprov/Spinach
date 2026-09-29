# examples/dnp_sol/cross_effect_field_scan_1.m

- MATLAB implementation: [examples/dnp_sol/cross_effect_field_scan_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/cross_effect_field_scan_1.m)

- **Call:** `cross_effect_field_scan_1()` (zero input arguments; the function plots the calculated curve and returns no explicit value).

## Purpose

Computes the steady-state proton magnetisation for a powder-averaged cross-effect DNP model under microwave irradiation, scanning applied-field offset. The source estimates hours of calculation time.

## Spin system and physical parameters

The source sets `sys.magnet=18.78` and uses `{'E','E','14N','1H'}`. The two electron g-tensors are `[2.0085 2.00605 2.00215]` (Euler angles `[pi/2 pi/3 pi/4]`) and `[2.00319 2.00319 2.00258]` (Euler angles `[pi/5 pi/6 pi/7]`). The `14N` quadrupolar entries are `[-1e6 -1e6 2e6]` with Euler angles `[0 0 0]`. Cartesian coordinates (explicitly labelled Angstrom) are `[0,0,0]`, `[12.80,0,0]`, unspecified for the `14N` (empty coordinate entry), and `[5.0,6.0,7.0]`. The electron-1/`14N` hyperfine entries are `[17.4e6 17.6e6 102e6]` with zero Euler angles. The electron-electron scalar coupling is assigned as `2*(-73e6)`.

The full `sphten-liouv` basis has no approximation. Relaxation uses `{'t1_t2'}` with `r1_rates={1e5,1e5,1e4,1e3}` and `r2_rates={1e7,1e7,1e5,1e4}`; `rlx_keep='diagonal'`, `equilibrium='zero'`, and `temperature=10` are specified. The source supplies no units for the relaxation rates, temperature, or the system-field assignment `18.78`.

## Calculation and observable

The ESR parameters set `mw_pwr=10e6`, `mw_frq=0`, electron irradiation spins `{'E'}`, `mw_oper=operator(spin_system,'Lx','E')/2`, and `ez_oper=operator(spin_system,'Lz','E')`. The script sets `parameters.fields=linspace(-0.08,0.04,256)`; this scanned offset axis is labelled Tesla in the plot. Powder averaging uses `rep_2ang_1600pts_sph`, method `backslash`, and `needs={'aniso_eq'}`. It evaluates `powder(spin_system,@dnp_field_scan,parameters,'esr')`, then plots `real(answer)` against the field offsets.

**Observable-label note:** the receiver is assigned as `state(spin_system,'Lz','1H')`, while the plot's vertical label says “$S_z$ expectation value on $^1$H”. The page reports both source statements without resolving their notation.

**Dependencies:** Spinach create/basis/state/operator, powder and plotting routines; the `dnp_field_scan` sequence and `rep_2ang_1600pts_sph` grid.
