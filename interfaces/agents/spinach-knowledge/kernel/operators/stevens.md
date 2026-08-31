# kernel/operators/stevens.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/stevens.m`
- Signature: `S=stevens(mult,k,q)`
- Total lines: 95

## Purpose

Extended Stevens operators. Syntax: S=stevens(mult,k,q)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(mult,k,q)`.
- Lines 33-34: Stockpile the ugly integer coefficients (there is no systematic formula); implemented by `C{13}=[1916006400 958003200 958003200 31933440 3991680 1995840 31680 15840 1584 264 24 12 1 ]`.
- Lines 48-49: Get Pauli matrices; implemented by `L=pauli(mult)`.
- Lines 51-52: Get the top state; implemented by `S=L.p^k`.
- Lines 54-55: Commute down; implemented by `for n=abs(q):(k-1)`.
- Lines 59-60: Distinguish odd and even indices; implemented by `if (~logical(mod(k,2)))&&(logical(mod(q,2)))`.
- Lines 66-67: Distinguish positive and negative indices; implemented by `if (q>0)||(q==0)`.

### Control flow inferred from the code

- Line 55: `for` loop over `n=abs(q):(k-1)`.
- Line 60: conditional branch on `(~logical(mod(k,2)))&&(logical(mod(q,2)))`.
- Line 67: conditional branch on `(q>0)||(q==0)`.

### Key state/data transformations

- Lines 34: computes `C{13}` using `C{13}=[1916006400 958003200 958003200 31933440 3991680 1995840 31680 15840 1584 264 24 12 1 ]`.
- Lines 35: computes `C{12}` using `C{12}=[319334400 79833600 79833600 13305600 2661120 23760 7920 1320 1320 22 22 1 ]`.
- Lines 36: computes `C{11}` using `C{11}=[14515200 7257600 1209600 604800 86400 2880 360 180 20 10 1 ]`.
- Lines 37: computes `C{10}` using `C{10}=[1451520 725700 725700 60480 60480 864 288 18 18 1 ]`.
- Lines 38: computes `C{9}` using `C{9} =[80640 40320 40320 6720 672 336 16 8 1 ]`.
- Lines 39: computes `C{8}` using `C{8} =[40320 5040 1680 168 168 14 14 1 ]`.
- Lines 40: computes `C{7}` using `C{7} =[2880 1440 360 60 12 6 1 ]`.
- Lines 41: computes `C{6}` using `C{6} =[480 240 240 10 10 1 ]`.
- Lines 42: computes `C{5}` using `C{5} =[48 24 8 4 1 ]`.
- Lines 43: computes `C{4}` using `C{4} =[24 6 6 1 ]`.
- Lines 44: computes `C{3}` using `C{3} =[4 2 1 ]`.
- Lines 45: computes `C{2}` using `C{2} =[2 1 ]`.
- Lines 46: computes `C{1}` using `C{1} =[1 ]`.
- Lines 49: computes `L` using `L=pauli(mult)`.
- Lines 52: computes `S` using `S=L.p^k`.
- Lines 61: computes `coeff` using `coeff=((-1)^(k-q))/C{k+1}(abs(q)+1)/2`.

### Local helper functions

- Line 76: `grumble()` — `function grumble(mult,k,q)`.
  - Representative operation: `if (~isnumeric(mult))||(~isreal(mult))|| (~isscalar(mult))||(mod(mult,1)~=0)||(mult<0)`.
  - Representative operation: `(~isscalar(mult))||(mod(mult,1)~=0)||(mult<0)`.

## Parameters / inputs

- mult -multiplicity of the spin in question
- k -Stevens operator rank, a non-negative
- integer
- q -Stevens operator projection, an
- integer between -k and k

## Outputs

- S -Stevens operator matrix
- Note: for historical reasons, the definition of Stevens operators
- is irregular and must rely on explicitly stockpiled coeffi-
- cients. Only ranks smaller or equal to 12 are available.

## Implementation structure

- Extended Stevens operators. Syntax:
- S=stevens(mult,k,q)
- mult -multiplicity of the spin in question
- k -Stevens operator rank, a non-negative
- integer
- q -Stevens operator projection, an
- integer between -k and k
- S -Stevens operator matrix
- Note: for historical reasons, the definition of Stevens operators
- is irregular and must rely on explicitly stockpiled coeffi-
- cients. Only ranks smaller or equal to 12 are available.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `logical()`, `isscalar()`.
