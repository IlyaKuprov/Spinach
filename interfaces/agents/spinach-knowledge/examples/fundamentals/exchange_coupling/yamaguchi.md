# examples/fundamentals/exchange_coupling/yamaguchi.m

- Signature: `yamaguchi()`

## Purpose

Estimates exchange coupling with the Yamaguchi equation from broken-symmetry DFT calculations for a bistrityl biradical with an alkynyl linker. The example is attributed to Olav Schiemann.

## Implementation

The script reads `biradical_singlet.log` and `biradical_triplet.log` with `gparse`, then calculates J = brokensymm(props_sing,props_trip). It displays J/1e9 under the label “Exchange coupling:” in GHz.
