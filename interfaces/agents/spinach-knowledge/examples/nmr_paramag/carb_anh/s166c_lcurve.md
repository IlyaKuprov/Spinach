# examples/nmr_paramag/carb_anh/s166c_lcurve.m

- Signature: `s166c_lcurve()`

## Purpose
L-curves for the S166C mutant dataset for human carbonic anhydrase II. The system and method are described in the [cited article](http://dx.doi.org/10.1039/c6sc03736d). A step-by-step tutorial is [available here](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).
## Physical / mathematical content


- The S166C example analyzes regularisation in the PCS inverse problem using an L-curve.

## Numerical / algorithmic content


- Sweeps 30 logarithmically spaced regularisation parameters, evaluates the PCS inverse problem in a parallel loop with GPU execution enabled, and uses L-curve analysis to suggest the regularisation parameter.

## Implementation structure


- Load experimental data
- Load susceptibility tensor
- Solver parameters
- Regularisation parameter array
- Result arrays
- Run a parallel loop
- L-curve analysis
