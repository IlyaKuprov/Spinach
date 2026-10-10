# examples/microfluidics/reacting_flow_nmr.m

## Coupled model

A spatial Diels–Alder reaction model on imported COMSOL mesh and velocity data couples advection-diffusion-reaction transport to proton NMR detection. The mesh crop is x=`[286.8, 287.5]` and y=`[576.0, 579.0]`, with listed mesh indices inactivated. The reaction records set the endo/exo channel rates to 2.0 and 1.0 L/(mol*s), respectively (the original numerical product assignments are preserved) and diffusion to `1e-7` without stating its unit. The initial field places 0.50 and 0.25 in reactant components at cells 1240 and 1246; their units are not given. Chemistry uses 501 steps of 20 seconds, while cellwise concentrations are interpolated with `makima` for spin evolution. Both stages use the kernel `kinetics` path. Removing all spins with `kill_spin` gives five concentration-only blocks for the first stage; `chem_concs` reads the propagated voxel populations. The frozen-rate stepping workflow is retained, but the kernel shares product unit arrival equally between reactants instead of assigning it entirely to one. Equal instantaneous mass-action derivatives therefore do not imply equal finite frozen steps: `test_cwdm_spatial` demonstrates this difference on a two-cell model. Full-chip history equivalence has not been established. The full spin-stage handle reads the prescribed concentrations from voxel unit coordinates and uses additive closure. Solvent retains three protons and T1/T2 relaxation but is unexcited and contributes no signal; its initial spatial population remains zero as in the original example.

The spin system uses `sys.magnet=14.1` (the file does not annotate a field unit) and a rectangular coil-selection mask defined by `287.0<x<287.3` and `577.0<y<577.5`. Proton acquisition parameters are offset 2328, sweep 3500, and 1024 steps; the file does not state units for offset or sweep, and sets the sampling interval to `1/sweep`. Unweighted `coil_state` vectors supply detection and longitudinal reference shapes; the propagated state also carries each voxel’s unit populations. It applies a `pi/2` excitation and advances the transport-, relaxation-, and concentration-dependent spin evolution with `step` using left/right interval generators.

## Observable

Acquisitions begin at every 25th chemistry-grid index (21 sampled starts, from 0 through 10,000 seconds). Each FID is apodised, zero-filled to 16384 points, Fourier transformed, and displayed as a real-intensity waterfall against time (seconds) and chemical shift (ppm); intensity is labelled a.u. The source header estimates days and says GPU execution is much faster; this is not a benchmark. The script defaults to `sys.enable={'greedy'}`; GPU transfers are conditional on `gpu` being enabled.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/reacting_flow_nmr.m

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
