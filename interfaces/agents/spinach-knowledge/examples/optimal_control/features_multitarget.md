# examples/optimal_control/features_multitarget.m

- Signature: `features_multitarget()`

## Purpose

Optimises a single pulse for two simultaneous state-transfer tasks: triplet–triplet to singlet–singlet (TT→SS) and triplet–singlet to singlet–triplet (TS→ST). The source evaluates the transfers over an ensemble of scalar couplings.

## Physical / mathematical content

The example uses proton and carbon controls, with x and y channels for each nucleus. It samples 11 J-coupling values spanning 13 to 16 and includes a state-norm (SNS) penalty in the optimisation objective.

## Numerical / algorithmic content

The source configures the lbfgs method and supplies grape_xy gradients to fmaxnewton. Both target transfers are assessed, rather than optimising a single target alone.

## Implementation structure

The script constructs the source and target states for TT→SS and TS→ST, sets the four proton/carbon x/y controls and 11-point coupling ensemble, and optimises with the SNS penalty. It then simulates the pulse and reports both transfer fidelities. The sampled J couplings span 13–16.
