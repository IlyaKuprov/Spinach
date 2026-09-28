# kernel/conventions/transforms/mhz2gauss.m

- Signature: `hfc_gauss=mhz2gauss(hfc_mhz,g)`

## Purpose

Converts hyperfine couplings from MHz to Gauss.

## Physical / mathematical content

With `muB=9.274009994e-24` and `hbar=1.054571628e-34`, the source uses `C=1e-10*g*muB/(hbar*2*pi)` and computes `hfc_gauss=hfc_mhz/C`.

## Numerical / algorithmic content

If omitted, `g` defaults to 2.0023193043622. The input coupling is divided by the frequency-per-field conversion factor C.

## Syntax

```matlab
hfc_gauss=mhz2gauss(hfc_mhz,g)
```

## Parameters / inputs

- `hfc_mhz` — real numeric array of hyperfine couplings in MHz.
- `g` — optional real numeric scalar electron g-factor; defaults to 2.0023193043622.

## Outputs

- `hfc_gauss` — converted hyperfine couplings in Gauss.

## Implementation structure

The function sets the default g-factor when needed, validates both inputs, then applies the conversion.
