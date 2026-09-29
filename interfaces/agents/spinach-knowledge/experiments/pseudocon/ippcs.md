# experiments/pseudocon/ippcs.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ippcs.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ippcs.m)

## Purpose

Fits a point-electron PCS model to measured pseudocontact shifts, estimating the paramagnetic-centre position and magnetic-susceptibility tensor. It is a nonlinear parameter-fitting utility, not a pulse-sequence simulator.

## Inputs and parameterisation

Call as `[mxyz,chi,pred_pcs,s_mxyz,s_chi] = ippcs(nxyz,mguess,expt_pcs)`. `nxyz` contains nuclear coordinates in Å as a real N-by-3 array (the implementation also accepts a cell array of such coordinate arrays). `mguess` is a real 1-by-3 initial estimate of the paramagnetic-centre coordinates in Å. `expt_pcs` is a real column of measured PCS values in ppm, one per nucleus.

The optimiser has eight free parameters: three centre coordinates and five independent susceptibility components in Å³. Those five values are assembled into a symmetric traceless tensor: the first row is `[p4 p5 p6]`, the second `[p5 p7 p8]`, and the third `[p6 p8 -p4-p7]`. The initial values for the five tensor parameters are each `0.1`. No parameter bounds or additional physical constraints are imposed by this fit.

## Fit and uncertainty estimates

The routine minimises the sum of squared residuals `expt_pcs - ppcs(nxyz,mxyz,chi)` using `fminunc`, central finite differences, and parallel function evaluations. The fitted PCS values are returned in `pred_pcs`. A good centre-position initial guess is important because the objective is nonlinear.

`jacobianest` estimates the residual Jacobian at the solution. The reported parameter standard deviations use the local least-squares covariance estimate `sdr^2 * inv(J'*J)`, where `sdr = sqrt(RSS/(N-8))`; the centre entries give `s_mxyz`. `s_chi` maps the five component estimates back onto a 3-by-3 layout, with the final diagonal estimate set to `sqrt(sp(4)^2+sp(7)^2)`. These are local linearised estimates from the fit, not independently validated confidence intervals.

## Outputs

`mxyz` is the fitted 1-by-3 centre coordinate, `chi` the fitted symmetric traceless susceptibility tensor in Å³, and `pred_pcs` the predicted ppm values at the input nuclei. `s_mxyz` and `s_chi` contain the corresponding coordinate and tensor-element standard-deviation estimates.

## References

- [Spin Dynamics Wiki: ippcs.m](https://spindynamics.org/wiki/index.php?title=ippcs.m)