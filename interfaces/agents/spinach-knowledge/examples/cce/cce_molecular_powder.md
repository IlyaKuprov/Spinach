# examples/cce/cce_molecular_powder.m

```matlab
[answer,diagnostics]=cce_molecular_powder();
```

Illustrative molecular Hahn echo of an electron and two deuterons. `clexpn` evaluates `@hahn_echo` at the three explicitly weighted orientations in `parameters.grid`, using independent cluster-local generators. The three axes illustrate tensor rotation, not powder convergence.

The three orientations are (0,0,0), (0,pi/2,0), (pi/2,pi/2,0), with equal weights. The 41-point time mesh spans 0 to 200 microseconds. This is a small illustration, not a converged powder grid or a published-material reproduction.

answer is a 41-by-1 normalised complex echo amplitude. diagnostics.reduction contains averaged order-resolved responses, validity, and absolute numerical bounds. diagnostics.orientations retains complete per-orientation evidence, including generalised Hahn raw observed densities and spectral diagnostics. Dense independent Kronecker/expm checks and the analytic zero-time limit reside in test_cce_examples.
