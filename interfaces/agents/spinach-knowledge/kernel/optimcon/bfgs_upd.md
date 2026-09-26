# kernel/optimcon/bfgs_upd.m

- Signature: `H=bfgs_upd(H,dx,dg)`

## Purpose

Performs a one-step BFGS Hessian update for maximisation using the argument and gradient increments from the previous step.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- This routine performs a single dense BFGS Hessian update for a maximisation problem. The sign conventions matter: the stored matrix approximates the negative Hessian so that solving H\g yields an ascent direction.
- The curvature test rejects bad secant pairs when dg^T dx does not have the sign and magnitude expected for locally concave behaviour. That protects the update from producing indefinite or numerically meaningless curvature models.

## Numerical / algorithmic content

- If no valid curvature information exists yet, the code scales an identity matrix using y^T y / y^T dx, a standard quasi-Newton initialisation that roughly matches curvature along the first accepted step.

## Syntax

```matlab
H=bfgs_upd(H,dx,dg)
```

## Parameters / inputs

- H -current BFGS approximation to the Hessian
- matrix corresponding to the *negative*
- Hessian of the objective; use [] on the
- first call
- dx -increment in arguments between the current
- and the previous step
- dg -increment in gradients between the current
- and the previous step

## Outputs

- H -updated BFGS approximation to the Hessian
- matrix corresponding to the *negative*
- Hessian of the objective

## Implementation structure

- Performs a one-step BFGS Hessian update for maximisation using
- the argument and gradient increments from the previous step.
- H=bfgs_upd(H,dx,dg)
- H -current BFGS approximation to the Hessian
- matrix corresponding to the *negative*
- Hessian of the objective; use [] on the
- first call
- dx -increment in arguments between the current
- and the previous step
- dg -increment in gradients between the current
- H -updated BFGS approximation to the Hessian
- Hessian of the objective
