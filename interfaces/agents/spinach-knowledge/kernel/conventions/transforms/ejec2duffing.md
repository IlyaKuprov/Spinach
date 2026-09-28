# kernel/conventions/transforms/ejec2duffing.m

- Signature: `[frq,anharm]=ejec2duffing(ej,ec)`

## Purpose

Converts the Josephson and charging energies of a transmon into the Duffing oscillator frequency and anharmonicity expected by the bosonic mode specification interface of create.m using the asymptotic transmon expressions (Koch et al., https://doi.org/ 10.1103/PhysRevA.76.042319): frq=sqrt(8*ej*ec)-ec, anharm=-ec

## Physical / mathematical content
- The function maps Josephson energy `ej` and charging energy `ec` to the Duffing transition frequency `sqrt(8*ej.*ec)-ec` and anharmonicity `-ec`.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

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
- Requires `ej` and `ec` to be real, finite, positive arrays of the same size.
- Computes `frq=sqrt(8*ej.*ec)-ec` and `anharm=-ec` elementwise.
- Warns when any `ej./ec` ratio is below 20.
