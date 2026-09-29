# examples/dnp_liq/odnp_liquid_1.m

- MATLAB implementation: [examples/dnp_liq/odnp_liquid_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/odnp_liquid_1.m)

- Signature: odnp_liquid_1()
- The source comment gives a calculation time of seconds.

## Purpose and run context

Demonstrates liquid-phase Overhauser DNP with continuous on-resonance CW electron irradiation; the source identifies Redfield theory as the treatment of dipolar cross-relaxation. Call odnp_liquid_1() with the Spinach liquid simulation routine and the dnp_time_dep callback available.

## Spin system and settings

There are two 1H spins and one electron, with sys.magnet=3.4. The proton Zeeman matrices are both diagonal with entries [5,5,5], and the electron matrix has diagonal entries [2.0023,2.0025,2.0027]. Coordinates are given in Å: proton 1 at [0,0,0], proton 2 at [0,2,0], and the electron at [0,0,1.5]. The complete basis uses sphten-liouv with no approximation.

Relaxation is {'redfield'}, equilibrium is dibari, the retained relaxation terms are secular, temperature is set to 298, and tau_c={10e-12}. The experiment requests electron spin E and rho_eq, and detects longitudinal Lz states for spins 1, 2, and 3. The electron microwave operator is Lx, its longitudinal operator is Lz, the offset is zero, and the microwave-power parameter is 2*pi*1e6. The time-step and step-count settings are 1e-6 and 1e3.

## Simulation output

The simulation call is liquid(spin_system,@dnp_time_dep,parameters,'esr'). The figure plots the electron longitudinal signal (answer row 3) and both proton longitudinal signals (rows 1–2) over 0–1000 microseconds. The proton legend labels the traces as 1.5 Å and 2.5 Å from the electron. The source plots these time traces; it does not provide a separate analytical enhancement value or a static-field comparison.
