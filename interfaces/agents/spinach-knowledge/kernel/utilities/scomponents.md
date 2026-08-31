# kernel/utilities/scomponents.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/scomponents.m`
- Signature: `sci=scomponents(A)`
- Total lines: 88

## Purpose

Strongly connected components of a graph, David Gleich's imple- mentation of Tarjan's algorithm:

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(A)`.
- Lines 31-32: Get the CSR indices; implemented by `[rp,ci]=sparse2csr(sparse(A))`.
- Lines 34-35: Run Tarjan's algorithm; implemented by `n=length(rp)-1; sci=zeros(n,1); cn=1`.

### Control flow inferred from the code

- Line 38: `for` loop over `sv=1:n`.
- Line 43: `while` loop over `rss>0`.
- Line 45: `while` loop over `ri<rp(v+1)`.
- Line 47: conditional branch on `root(w)==0`.
- Line 54: `for` loop over `ri=rp(v):rp(v+1)-1`.
- Line 56: conditional branch on `sci(w)==-1`.
- Line 57: conditional branch on `dt(root(v))>dt(root(w))`.
- Line 62: conditional branch on `root(v)==v`.
- Line 63: `while` loop over `css>0`.
- Line 65: conditional branch on `w==v, break; end`.

### Key state/data transformations

- Lines 32: computes `[rp,ci]` using `[rp,ci]=sparse2csr(sparse(A))`.
- Lines 35: computes `n` using `n=length(rp)-1; sci=zeros(n,1); cn=1`.
- Lines 36: computes `root` using `root=zeros(n,1); dt=zeros(n,1); t=0`.
- Lines 37: computes `cs` using `cs=zeros(n,1); css=0; rs=zeros(2*n,1); rss=0`.
- Lines 39: computes `v` using `v=sv; if root(v)>0, continue; end`.
- Lines 40: computes `rss` using `rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=rp(v)`.
- Lines 41: computes `root(v)` using `root(v)=v; sci(v)=-1; dt(v)=t; t=t+1`.
- Lines 42: computes `css` using `css=css+1; cs(css)=v`.
- Lines 46: computes `w` using `w=ci(ri); ri=ri+1`.
- Lines 48: computes `root(w)` using `root(w)=w; sci(w)=-1; dt(w)=t; t=t+1`.
- Lines 50: computes `rs(2*rss-1)` using `rs(2*rss-1)=v; rs(2*rss)=ri`.
- Lines 67: computes `cn` using `cn=cn+1`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(A)`. The Monks of Cool, whose tiny and exclusive monastery is hidden in a really cool and laid-back valley in the lower Ramtops, have a passing-
  - Representative operation: `if (~islogical(A))||(~ismatrix(A))||(size(A,1)~=size(A,2))`.
  - Representative operation: `error('the input must be a square logical matrix.')`.

## Syntax

```matlab
sci=scomponents(A)
```

## Parameters / inputs

- A -a logical square matrix with 1 for the
- connected nodes in the graph

## Outputs

- sci -a column vector with integers that spe-
- cify the strongly conected component
- that each node of the graph belongs to

## Implementation structure

- Strongly connected components of a graph, David Gleich's imple-
- mentation of Tarjan's algorithm:
- sci=scomponents(A)
- A - a logical square matrix with 1 for the
- connected nodes in the graph
- sci - a column vector with integers that spe-
- cify the strongly conected component
- that each node of the graph belongs to
- Check consistency
- Get the CSR indices
- Run Tarjan's algorithm
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sparse2csr()`, `root()`, `sci()`, `islogical()`, `ismatrix()`.
