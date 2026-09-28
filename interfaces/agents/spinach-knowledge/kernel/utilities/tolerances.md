# kernel/utilities/tolerances.m

- Signature: `[spin_system,sys]=tolerances(spin_system,sys)`

## Purpose

Sets numerical accuracy cut-offs and fundamental constants used by the Spinach kernel. The routine transfers tolerance settings from `sys.tols` into `spin_system.tols`; direct calls and modifications are discouraged, and accuracy settings should normally be supplied through `sys.tols` in the input preparation.

## Physical / mathematical content

This function configures numerical tolerances and constants; it does not construct interaction tensors or spin-dynamics propagators.

## Numerical / algorithmic content

After checking input consistency, the function applies user-specified tolerance fields, supplies defaults, handles the enabled paranoia options, reports the selected values, and rejects unrecognized tolerance options.

## Parameters / inputs

- `spin_system` — Spinach system description object.
- `sys` — system specification object described in the input preparation section of the manual.

## Outputs

- `spin_system` — updated system description object.
- `sys` — system specification structure with the tolerance substructure parsed out.

## Implementation structure

- Checks consistency with `grumble`.
- Reads configured tolerances from `sys.tols`, applies defaults and paranoia settings, and reports the resulting values.
- Source documentation: <https://spindynamics.org/wiki/index.php?title=tolerances.m>
