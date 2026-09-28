# examples/optimal_control/bloch_siegert/bloch_siegert_a.m

- Signature: `bloch_siegert_a()`

## Purpose

Bloch-Siegert shift compensation demo. It optimises a 90-degree pulse taking (L_z) to (L_x) for a single on-resonance (^{1}mathrm{H}) spin, then compares performance with and without Bloch-Siegert (BSS) correction as control power varies. Calculation time: minutes.

## Physical / mathematical content

- The source sets `sys.magnet=14.1`, one `1H` isotope, and zero scalar coupling. It builds normalized (L_z) initial and (L_x) target states. The example illustrates how an unaccounted Bloch-Siegert shift can reduce pulse fidelity at higher control power.

## Numerical / algorithmic content

- The optimizer uses L-BFGS (`max_iter=500`, `tol_x=1e-4`) and a 50-slice GRAPE-XY pulse; each slice duration is set to `(pi/100)/control.pwr_levels`. It sweeps 20 powers from (10^{-3}) to 1 times the proton Zeeman frequency, using the same random initial guess for each comparison. At each power, one pulse is designed with BSS disabled and one with BSS enabled; both are evaluated using the BSS-enabled setting, and terminal infidelity is plotted against relative control power.

## Implementation structure

- Create the one-spin system and sphten-liouv basis; build and normalize the initial and target states; obtain the (L_x)/(L_y) control operators and drift Hamiltonian; configure the optimizer; sweep power, optimize both BSS settings, evaluate both pulses under BSS, and plot the two infidelity curves.
