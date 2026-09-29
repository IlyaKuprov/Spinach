# kernel/integrity/exorcise.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/integrity/exorcise.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=exorcise.m).

- Signature: `exorcise(mode)`

## Purpose

Scans MATLAB files in the Spinach distribution for source-convention violations. It is a source-integrity check, not a numerical, chemical-flow, or line-shape routine.

## Scan and checks

The source recursively enumerates `.m` files under `kernel`, `interfaces`, `experiments`, and `etc`, excludes the listed foreign package directory `jsonlab-1.5`, and randomises the file order. It checks formatting (including consecutive blank lines, tabs, and the required file ending), function structure and the `grumble` consistency check, the length and sections of the introductory documentation, a Wiki link, portable `filesep` path construction, explicit norm types, and an `otherwise` branch in each `switch`. It also flags `disp` use when `spin_system` is available. Source markers `#NGRUM`, `#NHEAD`, `#NWIKI`, and `#NORMOK` mark the corresponding documented exceptions.

For source violations, the first failing file encountered is opened with MATLAB `edit` and the routine raises an error; because the scan order is randomised, this is not an exhaustive report of every failing file. The online mode additionally checks whether the URL in a file's header can be reached as a documentation page; offline mode skips that Wiki check.

## Input and output

`mode` must be a character string equal to `'online'` or `'offline'`; other values fail in `grumble(mode)`. The routine returns no value. It does not compute a physical quantity, so no equation, normalisation, output shape, or frequency units apply.

## Source guard

The mode/type validation is the public-input guard. Source-level gates are enforced as encountered and halt at the first detected violation after opening the file in the editor.

Related checks: [`existentials.m`](./existentials.md) checks startup prerequisites and path visibility; [`patrol.m`](./patrol.md) checks and executes selected examples.
