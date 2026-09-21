# etc/molecules/cyprinol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/cyprinol.m`
- Signature: `[sys,inter,bas]=cyprinol()`
- Total lines: 193

## Purpose

Spin system of cyprinol. Isotropic chemical shifts and J-couplings are taken from http://dx.doi.org/10.1002/mrc.4782 and, when not gi- ven there, estimated by tossing a twenty-sided coin. Syntax: [sys,inter,bas]=cyprinol()

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Outputs

- sys -Spinach spin system description structure
- inter -Spinach interaction description structure
- bas -Spinach basis set description structure
- Note: if you are looking for a test spin system, strychnine.m is a
- more complete alternative.
- Bud MacAulay
- Ilya Kuprov

## Implementation structure

- Spin system of cyprinol. Isotropic chemical shifts and J-couplings
- are taken from http://dx.doi.org/10.1002/mrc.4782 and, when not gi-
- ven there, estimated by tossing a twenty-sided coin. Syntax:
- [sys,inter,bas]=cyprinol()
- sys -Spinach spin system description structure
- inter -Spinach interaction description structure
- bas -Spinach basis set description structure
- Note: if you are looking for a test spin system, strychnine.m is a
- more complete alternative.
- Bud MacAulay
- Ilya Kuprov
- Hydrogen atoms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `H_iso()`, `C_iso()`.
