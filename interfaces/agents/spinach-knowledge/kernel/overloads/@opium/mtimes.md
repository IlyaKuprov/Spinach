# kernel/overloads/@opium/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@opium/mtimes.m`
- Signature: `c=mtimes(a,b)`
- Total lines: 85

## Purpose

Matrix products involving an OPIUM object. Syntax: c=mtimes(a,b)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: When A is a scalar; implemented by `if ~isa(a,'opium')&&isnumeric(a)&&isscalar(a)`.
- Lines 22-23: Return opium multiplied by A; implemented by `c=b; c.coeff=a*c.coeff; return`.
- Lines 27-28: When A is not a scalar; implemented by `if ~isa(a,'opium')&&isnumeric(a)`.
- Lines 30-31: Check dimension; implemented by `if size(a,2)~=size(b,1)`.
- Lines 35-36: Return A multiplied by opium; implemented by `c=b.coeff*a; return`.
- Lines 40-41: When B is a scalar; implemented by `if ~isa(b,'opium')&&isnumeric(b)&&isscalar(b)`.
- Lines 43-44: Return opium multiplied by A; implemented by `c=a; c.coeff=b*c.coeff; return`.
- Lines 48-49: When B is not a scalar; implemented by `if ~isa(b,'opium')&&isnumeric(b)`.
- Lines 56-57: Return A multiplied by opium; implemented by `c=a.coeff*b; return`.
- Lines 61-62: When both are opia; implemented by `if isa(a,'opium')&&isa(b,'opium')`.
- Lines 69-70: Update the coefficient; implemented by `c=a; c.coeff=a.coeff*b.coeff; return`.
- Lines 74-75: Complain and bomb out; implemented by `error('operands must be either numeric or opium objects.')`.

### Control flow inferred from the code

- Line 20: conditional branch on `~isa(a,'opium')&&isnumeric(a)&&isscalar(a)`.
- Line 28: conditional branch on `~isa(a,'opium')&&isnumeric(a)`.
- Line 31: conditional branch on `size(a,2)~=size(b,1)`.
- Line 41: conditional branch on `~isa(b,'opium')&&isnumeric(b)&&isscalar(b)`.
- Line 49: conditional branch on `~isa(b,'opium')&&isnumeric(b)`.
- Line 52: conditional branch on `size(a,2)~=size(b,1)`.
- Line 62: conditional branch on `isa(a,'opium')&&isa(b,'opium')`.
- Line 65: conditional branch on `size(a,2)~=size(b,1)`.

### Key state/data transformations

- Lines 23: computes `c` using `c=b; c.coeff=a*c.coeff; return`.

## Parameters / inputs

- a,b -opia or numerical arrays

## Outputs

- c -multiplication result

## Implementation structure

- Matrix products involving an OPIUM object. Syntax:
- c=mtimes(a,b)
- a,b -opia or numerical arrays
- c -multiplication result
- When A is a scalar
- Return opium multiplied by A
- When A is not a scalar
- Check dimension
- Return A multiplied by opium
- When B is a scalar
- When B is not a scalar
- When both are opia

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isscalar()`.
