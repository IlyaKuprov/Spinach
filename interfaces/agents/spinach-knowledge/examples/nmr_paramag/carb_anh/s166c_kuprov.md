# examples/nmr_paramag/carb_anh/s166c_kuprov.m

- Signature: `s166c_kuprov()`

## Purpose
Distributed fit for the S166C mutant dataset for human carbonic anhydrase II. The system and method are described in the [cited article](http://dx.doi.org/10.1039/c6sc03736d). A step-by-step tutorial is [available here](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).
## Physical / mathematical content


- The S166C example fits a distributed PCS source and derives an effective susceptibility tensor from that distribution.

## Numerical / algorithmic content

## Implementation structure


- Load experimental data
- Load susceptibility tensor
- Set inverse problem parameters
- Solve and refine the grid
- Get the new susceptibility tensor
