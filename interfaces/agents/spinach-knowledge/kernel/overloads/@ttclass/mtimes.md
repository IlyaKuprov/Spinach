# kernel/overloads/@ttclass/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/mtimes.m`
- Signature: `c=mtimes(a,b)`
- Total lines: 171

## Purpose

Performs tensor train multiplication followed by a shrink. Syntax: c=mtimes(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Decide the type combination; implemented by `if isa(a,'ttclass')&&isa(b,'double')&&isscalar(b)`.
- Lines 25-26: Multiply tensor train by a scalar from the right; implemented by `c=a; c.coeff=b*c.coeff; c.tolerance=abs(b)*c.tolerance`.
- Lines 30-31: Read sizes and ranks of the operands; implemented by `[a_ncores,a_ntrains]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a); b_sizes=size(b)`.
- Lines 33-34: Check consistency; implemented by `if (prod(a_sizes(:,2))~=b_sizes(1))`.
- Lines 38-39: Preallocate result; implemented by `mc=prod(a_sizes(:,1)); nc=b_sizes(2); c=zeros(mc,nc)`.
- Lines 41-42: Loop over the buffers of the operands; implemented by `for na=1:a_ntrains`.
- Lines 44-45: Set current vector as the right-hand side; implemented by `d=b`.
- Lines 47-48: Loop over the cores and perform multiplication; implemented by `for ka=a_ncores:(-1):1`.
- Lines 50-51: Reshape the current vector; implemented by `d=reshape(d,[a_ranks(ka+1,na)*a_sizes(ka,2),prod(a_sizes(1:ka-1,2))*prod(a_sizes(ka+1:a_ncores,1))*nc])`.
- Lines 53-54: Extract the core of the first operand and reshape it appropriately; implemented by `a_core=a.cores{ka,na}`.
- Lines 59-60: Perform the contraction; implemented by `d=a_core*d`.
- Lines 66-67: Accumulate the result; implemented by `c=c+a.coeff(na)*reshape(d,[mc,nc])`.
- Lines 72-73: Multiply tensor train by a scalar from the left; implemented by `c=b; c.coeff=a*c.coeff; c.tolerance=abs(a)*c.tolerance`.
- Lines 77-78: Read sizes and ranks of the operands; implemented by `[a_ncores,a_ntrains]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a)`.
- Lines 81-82: Check consistency; implemented by `if (a_ncores~=b_ncores)||(~all(a_sizes(:,2)==b_sizes(:,1)))`.
- Lines 86-87: Preallocate the result; implemented by `new_cores=cell(a_ncores,a_ntrains,b_ntrains)`.
- Lines 91-92: Loop over the buffers of the operands; implemented by `for nb=1:b_ntrains`.
- Lines 95-96: Loop over the cores; implemented by `for nc=1:a_ncores`.

### Control flow inferred from the code

- Line 23: conditional branch on `isa(a,'ttclass')&&isa(b,'double')&&isscalar(b)`.
- Line 34: conditional branch on `(prod(a_sizes(:,2))~=b_sizes(1))`.
- Line 42: `for` loop over `na=1:a_ntrains`.
- Line 48: `for` loop over `ka=a_ncores:(-1):1`.
- Line 82: conditional branch on `(a_ncores~=b_ncores)||(~all(a_sizes(:,2)==b_sizes(:,1)))`.
- Line 92: `for` loop over `nb=1:b_ntrains`.
- Line 93: `for` loop over `na=1:a_ntrains`.
- Line 96: `for` loop over `nc=1:a_ncores`.
- Line 115: conditional branch on `a.coeff(na)==0 || b.coeff(nb)==0`.
- Line 145: conditional branch on `~isempty(pos)`.

### Key state/data transformations

- Lines 26: computes `c` using `c=a; c.coeff=b*c.coeff; c.tolerance=abs(b)*c.tolerance`.
- Lines 31: computes `[a_ncores,a_ntrains]` using `[a_ncores,a_ntrains]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a); b_sizes=size(b)`.
- Lines 39: computes `mc` using `mc=prod(a_sizes(:,1)); nc=b_sizes(2); c=zeros(mc,nc)`.
- Lines 45: computes `d` using `d=b`.
- Lines 54: computes `a_core` using `a_core=a.cores{ka,na}`.
- Lines 79: computes `[b_ncores,b_ntrains]` using `[b_ncores,b_ntrains]=size(b.cores); b_ranks=ranks(b); b_sizes=sizes(b)`.
- Lines 87: computes `new_cores` using `new_cores=cell(a_ncores,a_ntrains,b_ntrains)`.
- Lines 88: computes `new_coeff` using `new_coeff=zeros(a_ntrains,b_ntrains)`.
- Lines 89: computes `new_tolerance` using `new_tolerance=zeros(a_ntrains,b_ntrains)`.
- Lines 99: computes `core_of_a` using `core_of_a=a.cores{nc,na}`.
- Lines 100: computes `core_of_b` using `core_of_b=b.cores{nc,nb}`.
- Lines 105: computes `core_of_c` using `core_of_c=reshape(core_of_a*core_of_b,[a_ranks(nc,na),a_ranks(nc+1,na),a_sizes(nc,1),b_sizes(nc,2),b_ranks(nc,nb),b_ranks(nc+1,nb)])`.
- Lines 109: computes `new_cores{nc,na,nb}` using `new_cores{nc,na,nb}=reshape(core_of_c,[a_ranks(nc,na)*b_ranks(nc,nb),a_sizes(nc,1),b_sizes(nc,2),a_ranks(nc+1,na)*b_ranks(nc+1,nb)])`.
- Lines 114: computes `new_coeff(na,nb)` using `new_coeff(na,nb)=a.coeff(na)*b.coeff(nb)`.
- Lines 118: computes `new_tolerance(na,nb)` using `new_tolerance(na,nb)=0`.
- Lines 123: computes `tt_current` using `tt_current=ttclass`.
- Lines 124: computes `tt_current.coeff` using `tt_current.coeff=a.coeff(na)*b.coeff(nb)`.
- Lines 125: computes `tt_current.cores` using `tt_current.cores=cell(a_ncores,1)`.

## Parameters / inputs

- a -a scalar or a tensor train
- b -a scalar, a tensor train, or a full matrix

## Outputs

- c -a tensor train object

## Implementation structure

- Performs tensor train multiplication followed by a shrink. Syntax:
- c=mtimes(a,b)
- a -a scalar or a tensor train
- b -a scalar, a tensor train, or a full matrix
- c -a tensor train object
- Decide the type combination
- Multiply tensor train by a scalar from the right
- Read sizes and ranks of the operands
- Check consistency
- Preallocate result
- Loop over the buffers of the operands
- Set current vector as the right-hand side

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isscalar()`, `ranks()`, `sizes()`, `a_sizes()`, `b_sizes()`, `a_ranks()`, `all()`, `b_ranks()`, `new_coeff()`, `new_tolerance()`, `new_cores()`, `unit_like()`, `shrink()`.
