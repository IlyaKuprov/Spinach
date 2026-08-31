# kernel/conventions/transforms/ejec2duffing.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/ejec2duffing.m`
- Signature: `[frq,anharm]=ejec2duffing(ej,ec)`
- Total lines: 72

## Purpose

Converts the Josephson and charging energies of a transmon into the Duffing oscillator frequency and anharmonicity expected by the bosonic mode specification interface of create.m using the asymptotic transmon expressions (Koch et al., https://doi.org/ 10.1103/PhysRevA.76.042319): frq=sqrt(8*ej*ec)-ec, anharm=-ec

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(ej,ec)`.
- Lines 44-45: Warn outside the transmon regime; implemented by `if any(ej./ec<20,'all')`.
- Lines 49-50: Run the conversion; implemented by `frq=sqrt(8*ej.*ec)-ec; anharm=-ec`.

### Control flow inferred from the code

- Line 45: conditional branch on `any(ej./ec<20,'all')`.

### Key state/data transformations

- Lines 50: computes `frq` using `frq=sqrt(8*ej.*ec)-ec; anharm=-ec`.

### Local helper functions

- Line 55: `grumble()` — `function grumble(ej,ec)`.
  - Representative operation: `if (~isnumeric(ej))||(~isreal(ej))||any(~isfinite(ej),'all')||any(ej<=0,'all')`.
  - Representative operation: `error('ej must be an array of positive real numbers.')`.

## Syntax

```matlab
[frq,anharm]=ejec2duffing(ej,ec)
```

## Parameters / inputs

- ej -Josephson energies in Hz (energy over the
- Planck constant), an array of positive
- real numbers
- ec -charging energies in Hz (energy over the
- Planck constant), an array of positive
- real numbers of the same size as ej

## Outputs

- frq -transition frequencies in Hz, to be placed
- into inter.modes.frqs
- anharm -Duffing anharmonicities in Hz, to be placed
- into inter.modes.anharms
- Note: the asymptotic expressions are only accurate deep in the
- transmon regime ej/ec>>1; a warning is issued when the
- ratio is smaller than 20.

## Implementation structure

- Converts the Josephson and charging energies of a transmon into
- the Duffing oscillator frequency and anharmonicity expected by
- the bosonic mode specification interface of create.m using the
- asymptotic transmon expressions (Koch et al., https://doi.org/
- 10.1103/PhysRevA.76.042319):
- frq=sqrt(8*ej*ec)-ec, anharm=-ec
- [frq,anharm]=ejec2duffing(ej,ec)
- ej -Josephson energies in Hz (energy over the
- Planck constant), an array of positive
- real numbers
- ec -charging energies in Hz (energy over the
- real numbers of the same size as ej

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `isequal()`.
