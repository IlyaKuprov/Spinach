# experiments/fieldsweep.m

- Signature: [spec,parameters]=fieldsweep(spin_system,parameters)

## Purpose and physical scope

Computes a simulated field-swept powder EPR spectrum by locating spin transitions at a fixed microwave frequency and integrating their orientation-dependent contributions over a spherical powder grid. The transition Hamiltonians come from the supplied Spinach system; electron-nuclear hyperfine effects appear only if represented in that system. This routine is not an ESEEM or ENDOR pulse sequence and does not calculate a measured spectrum.

## Inputs and parameters

- parameters.grid: initial spherical grid name. The source recommends a non-symmetric grid to avoid transition degeneracies and gives rep_2ang_100pts_sph as a starting example.
- parameters.spins: one-element cell array naming the spin coupled to the microwave field; the source example is {'E'}.
- parameters.mw_freq: microwave frequency in Hz.
- parameters.fwhm: Lorentzian line full width at half maximum in tesla.
- parameters.window: two-element field interval [Bmin Bmax] in tesla.
- parameters.npoints: number of field samples.
- parameters.tm_tol: relative transition-moment tolerance; 0.01 is the source's suggested starting value.
- parameters.rspt_order: perturbation-theory order for eigenfields; 2 is suggested, while Inf requests exact diagonalisation.
- parameters.int_tol: powder-integration tolerance; the source gives no numeric setting.
- Spinach system field: set the system magnetic field to 1 tesla, irrespective of the sweep window, as required by the source.

## Calculation and returned axes

The function obtains coupling and Zeeman Hamiltonian terms, constructs the unweighted microwave Lx operator vector with `coil_state` for the selected spin, loads the named spherical grid, and forms its convex hull. It builds parameters.b_axis as linspace(window(1),window(2),npoints). At each grid vertex it calls eigenfields for the specified microwave frequency and orientation; it then integrates contributions triangle by triangle with the recursive Voitlander integrator and sums the triangle spectra. The source describes the line shape as Lorentzian.

The first output spec is the field-sampled spectrum, with its field coordinates in the returned parameters.b_axis (tesla). The second output is the updated parameter structure, including that axis; the source signature does not return b_axis as a separate first or second output.

## Scope and limitations

This is a computational powder integration over a finite spherical grid, with a potentially expensive eigenfields calculation. The source's example settings are starting points, not validated accuracy guarantees; increase grid density or adjust tolerances as required by convergence. Use rspt_order=Inf when exact diagonalisation is desired instead of the suggested perturbative order 2. The function must be called directly, without an experiment context.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/fieldsweep.m
