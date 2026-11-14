# kernel/cce/cce_validate.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_validate.m)

## Purpose

Shared public CCE optional-field contract used by preparation, evaluators and parameter-consuming helpers.

## Behaviour

A finite `meanfield_radius` requires finite real central and exterior coordinates; missing geometry raises `cce:coordinates` before distance filtering. The radius is nonnegative, and `Inf` retains the unbounded default without requiring coordinates. Included cluster sites need no coordinates for this cutoff, and an empty exterior requires no distance calculation.

Rejects unknown top-level and nested fields, reserved extensions, inactive physical options, incompatible bath-state and graph options, and unsupported propagation preferences. `parameters.controls` is rejected even when empty: pulse events belong in `sequence`, and bandwidth selection uses `deer`. Helper-specific dimensions and values remain checked by each helper.

## Inputs and outputs

`cce_validate(parameters)`; evaluators additionally supply `(parameters,model,sites,context)` for geometry checks accepts one scalar structure and returns no value. A failed contract raises an informative error naming the unsupported field.

The actual second-order projection API remains top-level `second_order`; it is restricted to conditional CCE/ME-CCE. Numerics accepts the implemented `eig` preference for closed/noise propagation and `stock` for ME propagation. It does not accept the former unused `numerics.reduction` alias.
