# examples/fundamentals/paramag_test.m

- Signature: `paramag_test()`

## Purpose

Checks two internal consistency relationships in the paramagnetic chemical-shift functions: the isotropic part of `xyz2pms` versus `ppcs`, and the tensor conversion through `xyz2hfc` and `hfc2pms` versus direct `xyz2pms` evaluation.

## Physical and mathematical content

The function generates a random symmetric 3-by-3 susceptibility matrix `chi` and random 3D coordinate vectors. In the first comparison it checks `ppcs(nxyz,sxyz,chi)` against `trace(xyz2pms(nxyz,sxyz,chi))/3`, with an absolute scalar difference below `1e-6`. In the second, it uses isotope `15N`, obtains a hyperfine tensor from `xyz2hfc(mxyz,nxyz,isotope)`, converts it with `hfc2pms(A,chi,isotope)`, and compares the returned paramagnetic-shift tensor with `xyz2pms(nxyz,mxyz,chi)`. The Frobenius-norm difference must be below `1e-6`.

## Callable context and assumptions

Call the zero-input MATLAB function `paramag_test()` from a Spinach checkout with the project functions on the MATLAB path. It returns no values and prints `Test 1 passed.` or `Test 2 passed.` only on the corresponding successful branch; otherwise it raises an error. The susceptibility tensor is made symmetric as `(chi+chi')/2`, but the source does not impose positive definiteness or set a random seed. The coordinates and susceptibility entries are random, so this is a conversion-consistency check over sampled inputs, not a fixed benchmark or evidence for a particular experimental system.

## Numerical content and limits

The source defines two thresholds of `1e-6` and no fixed numerical outputs. This entry describes the assertions encoded in the function and does not state that they have been run or passed.

## Source

[`examples/fundamentals/paramag_test.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/paramag_test.m)
