# kernel/overloads/@ttclass/kron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/kron.m`
- Signature: `c=kron(a,b)`
- Total lines: 67

## Purpose

Kronecker product of two matrices in a tensor train format. Syntax: c=kron(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Shrink a and b before going any further; implemented by `a=shrink(a); b=shrink(b)`.
- Lines 29-30: Read sizes and ranks of the operands; implemented by `[a_ncores,~]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a)`.
- Lines 33-34: Check consistency; implemented by `if a_ncores~=b_ncores`.
- Lines 38-39: Preallocate the result; implemented by `new_cores=cell(a_ncores,1)`.
- Lines 41-42: Kron the cores; implemented by `for nc=1:a_ncores`.
- Lines 49-50: Write the output data structure; implemented by `c=ttclass; c.coeff=b.coeff*a.coeff`.
- Lines 53-54: The result is exactly zero; implemented by `c.tolerance=0`.
- Lines 56-57: The relative tolerances sum up; implemented by `c.tolerance=abs(a.coeff)*b.tolerance+abs(b.coeff)*a.tolerance`.

### Control flow inferred from the code

- Line 34: conditional branch on `a_ncores~=b_ncores`.
- Line 42: `for` loop over `nc=1:a_ncores`.
- Line 52: conditional branch on `a.coeff==0 || b.coeff==0`.

### Key state/data transformations

- Lines 27: computes `a` using `a=shrink(a); b=shrink(b)`.
- Lines 30: computes `[a_ncores,~]` using `[a_ncores,~]=size(a.cores); a_ranks=ranks(a); a_sizes=sizes(a)`.
- Lines 31: computes `[b_ncores,~]` using `[b_ncores,~]=size(b.cores); b_ranks=ranks(b); b_sizes=sizes(b)`.
- Lines 39: computes `new_cores` using `new_cores=cell(a_ncores,1)`.
- Lines 43: computes `core_of_a` using `core_of_a=reshape(a.cores{nc,1},[a_ranks(nc,1)*a_sizes(nc,1),a_sizes(nc,2)*a_ranks(nc+1,1)])`.
- Lines 44: computes `core_of_b` using `core_of_b=reshape(b.cores{nc,1},[b_ranks(nc,1)*b_sizes(nc,1),b_sizes(nc,2)*b_ranks(nc+1,1)])`.
- Lines 45: computes `core_of_c` using `core_of_c=reshape(kron(core_of_b,core_of_a),[a_ranks(nc,1),a_sizes(nc,1),b_ranks(nc,1),b_sizes(nc,1),a_sizes(nc,2),a_ranks(nc+1,1),b_sizes(nc,2),b_ranks(nc+1,1)])`.
- Lines 46: computes `new_cores{nc,1}` using `new_cores{nc,1}=reshape(permute(core_of_c,[1 3 4 2 7 5 6 8]),[a_ranks(nc,1)*b_ranks(nc,1),b_sizes(nc,1)*a_sizes(nc,1),b_sizes(nc,2)*a_sizes(nc,2),a_ranks(nc+1,1)*b_ranks…`.
- Lines 50: computes `c` using `c=ttclass; c.coeff=b.coeff*a.coeff`.
- Lines 51: computes `c.cores` using `c.cores=reshape(new_cores,[a_ncores,1])`.
- Lines 54: computes `c.tolerance` using `c.tolerance=0`.

## Parameters / inputs

- a,b -tensor train objects

## Outputs

- c -a tensor trin object
- WARNING: the result is not the same as the flat matrix Kronecker pro-
- duct (it is a row and column permutation away from it), but
- the resulting order of elements is consistent with the out-
- put of the tensor train vectorization (ttclass/vec) operati-
- on output.

## Implementation structure

- Kronecker product of two matrices in a tensor train format. Syntax:
- c=kron(a,b)
- a,b -tensor train objects
- c -a tensor trin object
- WARNING: the result is not the same as the flat matrix Kronecker pro-
- duct (it is a row and column permutation away from it), but
- the resulting order of elements is consistent with the out-
- put of the tensor train vectorization (ttclass/vec) operati-
- on output.
- Shrink a and b before going any further
- Read sizes and ranks of the operands
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `shrink()`, `ranks()`, `sizes()`, `a_ranks()`, `a_sizes()`, `b_ranks()`, `b_sizes()`.
