# tests/kernel/test_slowpass_sat.m

Regression checks for saturated continuous-wave detection. Compares both complex quadratures with analytic one-electron Bloch solutions from zero drive to strong saturation, with a nonzero resonance offset, and with a native propagator fixed point. Checks identity-containing detectors, ignored rho0, independently conserved non-unit concentrations including zero population, liquid/gridfree isotropic agreement, anisotropic rotational diffusion, and the absence of a ZTE call.

A sparse space-times-spin fixture above 100000 active coordinates exercises the production GMRES branch against the spin-only solution. Negative cases cover invalid nutation frequencies, nonfinite sweep endpoints, unthermalised relaxation, and population loss. This is not a substitute for saturated nitroxide spatial-basis or field-grid convergence.
