# Zeeman chemical maps and Hilbert matrix derivatives

`test_cwdm_zeeman_kin` compares production maps against independent physical matrix formulas. A reversed-spin A+B-to-C product retains complex coherences; its additive counterpart retains the identity once and drops cross-reactant polarisation products. The tests cover zero reactant population, separate spatial voxels, time rates, RKMK4 propagation against an independent product-state ODE, partial trace from a correlated two-spin singlet, and unpolarised arrival into an unmatched product spin.

Hilbert first-order exchange is compared with explicit matrix drains/fills and integrated against a full-component matrix exponential. `chem_concs` is checked against local traces. Hilbert IME is tested for stationarity, matrix action, trace preservation, and weighted-target rejection. Cross-species matrix coherences must be rejected; time-dependent first-order matrix rates are exercised. Named and user-superoperator singlet selectors must give the same Haberkorn loss and spin-free product arrival in both density-matrix Zeeman formalisms. Hilbert mass-action rejections assert the exact identifier and complete message for both closures.

Algebraic whole-array comparisons use 1e-12 absolute tolerance; integrated exchange uses 1e-10. Expected maps are not generated from spherical-tensor descriptors.
