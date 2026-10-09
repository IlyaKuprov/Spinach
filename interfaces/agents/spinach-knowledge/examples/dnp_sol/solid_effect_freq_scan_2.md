# examples/dnp_sol/solid_effect_freq_scan_2.m

- Signature: `solid_effect_freq_scan_2()`
- Source: [`examples/dnp_sol/solid_effect_freq_scan_2.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/solid_effect_freq_scan_2.m)

## Purpose

Scans the microwave frequency in a powder-averaged, laboratory-frame steady-state DNP calculation for a single `15N-labelled` urea spin system coupled to one electron. The example compares two widely separated frequency windows and observes longitudinal polarisation on both proton and nitrogen channels. The source estimates a calculation time of hours.

## Spin system and model

The source sets `sys.magnet=3.4` (the source calls this a magnetic field but does not state its unit) and orders the seven spins as `E, 15N, 1H, 1H, 15N, 1H, 1H`. It assigns the following literal coordinate vectors, in that order: electron [0, 0, 10.14358975]; nitrogen [-0.07640311, 1.16112702, -0.61556225]; proton [0.08533754, 1.99241453, -0.06489225]; proton [0.38824423, 1.16155815, -1.51333625]; nitrogen [0.07640311, -1.16112702, -0.61556225]; proton [-0.08533754, -1.99241453, -0.06489225]; proton [-0.38824423, -1.16155815, -1.51333625]. The source describes one labelled urea molecule at a selected orientation and electron distance, but does not annotate coordinate units.

The basis uses `sphten-liouv`, `IK-0`, interaction level 4, and projections [-2, -1, 0, +1, +2]. Relaxation is the Weizmann model with secular retention and zero equilibrium. The configured rates are `weiz_r1e=1e2`, `weiz_r1n=0.1`, `weiz_r2e=1e5`, and `weiz_r2n=1e3`; both 7-by-7 dipolar-rate arrays are set to `1e-3*ones(7,7)`. The source sets temperature to 4.2; it does not annotate units for these rates or temperature.

## Scan and computation

The function builds the spin system and basis with `create` and `basis`, then calls `powder(spin_system,@dnp_freq_scan,parameters,'esr')`. The electron channel, microwave power `2*pi*100e3`, two longitudinal detection operators (1H and 15N), the electron microwave and Zeeman operators, and `g_ref=spin_system.tols.freeg` are supplied in the parameters. The 200-point frequency vector concatenates 100 points from 144.0 to 145.5 MHz and 100 points from 14.0 to 15.5 MHz. Powder averaging uses `rep_2ang_100pts_sph`, method `lvn-backs`, and the `aniso_eq` requirement. The callback is Spinach's [`dnp_freq_scan`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/dnp_freq_scan.m).

## Result and scope

The returned array is plotted in four panels: 1H then 15N longitudinal expectation values over each of the two frequency windows, using the real parts of the corresponding result column and frequency block. The axes are labelled microwave frequency in MHz; no output array or figure is saved by the function. The source supplies no computed polarisation values or performance result beyond its hours estimate.
