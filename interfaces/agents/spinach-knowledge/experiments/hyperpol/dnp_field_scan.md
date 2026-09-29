# experiments/hyperpol/dnp_field_scan.m

- MATLAB implementation: [experiments/hyperpol/dnp_field_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/dnp_field_scan.m)

- Signature: `dnp=dnp_field_scan(spin_system,parameters,H,R,K)`

## Purpose and physical scope

This routine solves a steady-state DNP problem at each supplied magnetic-field offset. It is a steady-state linear-system calculation, not time propagation or a pulse-sequence controller. Microwave driving and the electron Zeeman offset act on the spin model supplied by the caller; the caller supplies the Hamiltonian, relaxation and kinetics superoperators. Electron-nuclear hyperfine interactions contribute only if they are present in that model. This is not an ESEEM/ENDOR sequence or an image-reconstruction routine.

## Inputs

- `parameters.mw_pwr`: microwave power in Hz; the implementation converts it to angular units in the Liouvillian.
- `parameters.mw_frq`: microwave-frequency offset in Hz from the free-electron frequency at the reference B0.
- `parameters.fields`: real vector of magnetic-field offsets from reference B0, in tesla.
- `parameters.rho0`: equilibrium state at reference B0.
- `parameters.coil`: one detection-state vector or a horizontal stack of detection states.
- `parameters.mw_oper`: microwave irradiation operator; `parameters.ez_oper`: electron `Lz` operator.
- `parameters.method`: `'backslash'` or `'gmres'`.
- H, R and K: Hamiltonian, relaxation and kinetics matrices supplied by the context function.

## Calculation and output axes

The routine forms the generator from H + `1i`*R + `1i`*K, adds the microwave and reference-frequency terms, and forms the source b = R*`rho0`. For each field it adds the electron Zeeman offset, solves the steady-state linear system, and projects the solution onto each detection state. The returned dnp has shape [numel(`parameters.fields`), size(`parameters.coil`,2)]: rows follow the supplied field vector and columns follow the `coil`-state stack. Values are expectation values and may be complex; the example plots their real part.

## Model limits

The source explicitly requires an unthermalised relaxation superoperator. It also assumes the equilibrium state and relaxation superoperator are the same at every field in the sweep, and warns against broad field sweeps. Use only where that fixed-reference assumption is suitable. Supported formalisms are `sphten-liouv` and `zeeman-liouv`.

## Source-coded numerical example

`examples/dnp_sol/solid_effect_field_scan_1.m` sets `parameters.mw_pwr`=1e5 Hz, `parameters.mw_frq`=-14e8 Hz, and `parameters.fields`=linspace(-0.08,+0.08,512) tesla. It detects the 15N `Lz` state in a gadolinium-containing solid-effect DNP example. These are simulation inputs, not measured values or a reported calculation result.

## Source and attribution

- Source: `experiments/hyperpol/dnp_field_scan.m`
- <https://spindynamics.org/wiki/index.php?title=dnp_field_scan.m>
- Source attribution: ilya.kuprov@weizmann.ac.il
