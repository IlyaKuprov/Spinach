# CWDM sequence probes

`test_cwdm_sequences` exercises the seven reviewed sequence paths. It checks finite zero-population relaxation rates, concentration-independent normalised ENDOR and microwave operators, linear concentration scaling of 2D/3D PRESS phantoms, and both 3D imaging acquisitions on homogeneous grids. The imaging input density and receiver are held fixed when stored concentrations change. This isolates the geometry contract; it is not a powder or spatial convergence test.

Unweighted operator-vector calls explicitly supply the `exact` method required by the four-argument `coil_state` API.
