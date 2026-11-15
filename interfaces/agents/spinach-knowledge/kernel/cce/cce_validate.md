# kernel/cce/cce_validate.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_validate.m)

## Purpose

Shared public CCE optional-field contract used by preparation, evaluators and parameter-consuming helpers.

## Behaviour

A finite `meanfield_radius` requires finite real central and exterior coordinates; missing geometry raises `cce:coordinates` before distance filtering. The radius is nonnegative, and `Inf` retains the unbounded default without requiring coordinates. Included cluster sites need no coordinates for this cutoff, and an empty exterior requires no distance calculation.

Rejects unknown top-level and nested fields, reserved extensions, inactive physical options, incompatible bath-state and graph options, and unsupported propagation preferences. `parameters.controls` is rejected even when empty: pulse events belong in `sequence`, and bandwidth selection uses `deer`. Helper-specific dimensions and values remain checked by each helper.

## Inputs and outputs

`cce_validate(parameters)`; evaluators additionally supply `(parameters,model,sites,context)` for geometry and state-context checks. The helper accepts a scalar parameter structure and returns no value. A failed contract raises an informative error naming the unsupported field.

The actual second-order projection API remains top-level `second_order`; it is restricted to conditional CCE/ME-CCE. Numerics accepts the implemented `eig` preference for closed/noise propagation and `stock` for ME propagation. It does not accept the former unused `numerics.reduction` alias.

Evaluator state contexts accept the canonical `cce_states` fields, the interlaced assignment and optional `id` identity metadata. Unsupported initial-state, projection and reference aliases are rejected with `cce:context`: use `context.central_rho` and `context.bath_rho` for density overrides, and `parameters.central.alpha/beta` for branch projections. The unused `parameters.initial_rho` alias is rejected on helper and driver paths alike.
