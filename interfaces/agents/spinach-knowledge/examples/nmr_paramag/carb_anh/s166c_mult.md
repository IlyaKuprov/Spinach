# examples/nmr_paramag/carb_anh/s166c_mult.m

- Signature: `s166c_mult()`

## Purpose
Multipolar fit for the S166C mutant dataset for human carbonic anhydrase II. The system and method are described in the [cited article](http://dx.doi.org/10.1039/c6sc03736d). A step-by-step tutorial is [available here](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).
## Physical / mathematical content


- The S166C example fits PCS data with a multipolar model and reports the susceptibility tensor and multipole-centre location.

## Numerical / algorithmic content

## Implementation structure


- Load experimental data
- Solve the inverse problem
- Plot experimental vs predicted PCS
- Report and save the parameters
