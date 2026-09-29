# examples/dnp_liq/odnp_liquid_3.m

- MATLAB implementation: [examples/dnp_liq/odnp_liquid_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/odnp_liquid_3.m)

- Signature: `odnp_liquid_3()`
- Calculation time: minutes
- [MATLAB source](../../../../../examples/dnp_liq/odnp_liquid_3.m)

## Purpose

Maps the steady-state proton longitudinal signal against microwave-frequency offset and magnetic field for a liquid-phase electron–nucleus DNP model. The example is constructed to show a high-field g–hyperfine cross-correlation effect. Rather than propagating a transient, it sets the time derivative of the inhomogeneous master equation to zero and solves for the steady-state density matrix.

## Spin system and relaxation model

The spins are one proton and one electron. The proton Zeeman eigenvalues are [15, 5, −20] ppm; the electron values are [2.00210, 2.00250, 2.00290] as dimensionless g factors. Their Euler angles are [0, 0, 0] and [pi/3, pi/4, pi/5], respectively. The isotropic hyperfine coupling is 20e6 Hz. The coordinates are (0, 0, 0) and (0, 0, 3.0) Angstrom for the proton and electron.

The basis is `sphten-liouv` with no approximation. Redfield relaxation uses zero equilibrium (required by this steady-state calculation), secular retention, temperature 298, and a 10 ps correlation time (identified in the source as TEMPOL in water). The relaxation-integration tolerance is set to `1e-10`, which the source marks as necessary for this calculation.

## Field and microwave scan

The field grid is `linspace(1,10,64)` Tesla. At each field, a `parfor` iteration sets the local system magnet, creates the spin system and basis, and constructs the proton `Lz` coil and electron `Lx` microwave and `Lz` offset operators. The ESR-context calculation is `liquid(spin_system,@dnp_freq_scan,locpar,'esr')`. Its parameters use the electron spin, `method='lvn-backs'`, `needs={'rho_eq'}`, and `g_ref=mean(inter.zeeman.eigs{2})`; microwave power is `2*pi*500e3`. The offset vector is `2*pi*linspace(-15,15,512)*1e6`, spanning −15 to +15 MHz in the plot's frequency-offset units.

## Output and scope

The result is a 512-by-64 array. The image plot displays its real part versus field and microwave-frequency offset; the colour bar is the steady-state proton longitudinal signal, labelled as `<H_Z>`. This is a steady-state scan over the specified single-spin-pair model and field/frequency grids; the source does not present a time-domain trajectory.
