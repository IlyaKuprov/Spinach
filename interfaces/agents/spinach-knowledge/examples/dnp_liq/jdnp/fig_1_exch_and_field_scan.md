# examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m

- Signature: `fig_1_exch_and_field_scan()`
- Reference: [Physical Chemistry Chemical Physics, DOI: 10.1039/D1CP04186J](https://doi.org/10.1039/d1cp04186j)
- Calculation time: hours (256 field values, each with 256 exchange-coupling simulations)

## Purpose

Maps the proton DNP enhancement after a 20 ms microwave pulse as a function of static field and electron–electron exchange coupling, for the spin system returned by `system_specification()`.

## Calculation

The script scans 256 fields from 0.25 to 3.0 T and 256 exchange couplings from -100 to +100 GHz. At each field it sets the microwave offset to the difference between the trityl and free-electron frequencies; for each exchange value it builds the Spinach system and basis, constructs the ESR Hamiltonian and relaxation superoperator, adds the microwave drive and offset, and propagates the thermal-equilibrium state for 20 ms. The plotted quantity is the detected proton magnetisation divided by its equilibrium value. The exchange loop uses `parfor`, and the field-by-exchange map is refreshed after each field row.
