# examples/dnp_liq/odnp_liquid_3.m

- Signature: `odnp_liquid_3()`
- Calculation time: minutes

## Purpose

Computes the steady-state proton longitudinal signal over microwave-frequency offset and magnetic field for a liquid-state electron–nucleus DNP model. The source describes the example as showing a high-field g–hyperfine cross-correlation effect.

## Spin system and steady-state calculation

The model contains one proton and one electron. It specifies anisotropic Zeeman tensors, a 20 MHz isotropic hyperfine coupling, and a 3 Å interspin separation for the rank-2 anisotropic dipolar interaction. The calculation uses a complete sphten-liouv basis, Redfield relaxation, zero equilibrium (as required by this steady-state setup), secular relaxation retention, 298 K, and a 10 ps correlation time. The source tightens the relaxation-integration tolerance to 1e-10.

The ESR-context `liquid` calculation calls `dnp_freq_scan` with method `lvn-backs` and 500 kHz microwave power. It scans 512 offsets from -15 to +15 MHz and 64 fields from 1 to 10 T; a `parfor` loop processes the field values. The output is a colour map of the real steady-state proton `Lz` signal versus field and microwave offset.
