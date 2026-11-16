# examples/cce/cce_molecular_powder.m

```matlab
[answer,diagnostics]=cce_molecular_powder();
```

Illustrative molecular Hahn echo of an electron and two deuterons. The example creates the physical system, declares a generalized second-order clexpn calculation and transverse electron preparation, binds the named hahn_echo experiment with context_task, and calls powder. No full-bath basis, explicit state or detector, pulse-event array, transition vector or manual orientation loop is assembled.

The unchanged three orientations are (0,0,0), (0,pi/2,0), (pi/2,pi/2,0), with equal weights. The unchanged 41-point time mesh spans 0 to 200 microseconds. This is a small illustration, not a converged powder grid or a published-material reproduction.

answer is a 41-by-1 normalized complex echo amplitude. diagnostics.reduction contains averaged order-resolved responses, validity and absolute numerical bounds. diagnostics.orientations retains complete per-orientation evidence, including generalized Hahn raw observed densities and spectral diagnostics. Dense independent Kronecker/expm checks and the analytic zero-time limit reside in test_cce_examples.
