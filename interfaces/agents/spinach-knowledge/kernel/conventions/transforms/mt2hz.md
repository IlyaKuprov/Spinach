# kernel/conventions/transforms/mt2hz.m

- Signature: `hfc_hz=mt2hz(hfc_mt,g)`

## Purpose

Converts hyperfine couplings from milliTesla (mT) to linear frequency in hertz (Hz).

## Physical / mathematical content

With `muB=9.274009994e-24` and `hbar=1.054571628e-34`, the conversion is `hfc_hz=hfc_mt*(1e-3*g*muB/(hbar*2*pi))`.

## Numerical / algorithmic content

If omitted, `g` defaults to 2.0023193043622. The conversion scales each input value by the same factor.

## Syntax

```matlab
hfc_hz=mt2hz(hfc_mt,g)
```

## Parameters / inputs

- `hfc_mt` — real numeric array of hyperfine couplings in mT.
- `g` — optional real numeric scalar electron g-factor; defaults to 2.0023193043622.

## Outputs

- `hfc_hz` — converted hyperfine couplings in Hz.

## Implementation structure

The function applies the default g-factor when it is omitted, validates the inputs, then multiplies by the conversion factor.
