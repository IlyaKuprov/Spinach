# examples/fundamentals/paramag_test.m

- Signature: `paramag_test()`

## Purpose

Check two equivalent calculations of paramagnetic chemical shifts: the isotropic shift from ppcs versus the trace of xyz2pms, and the tensor from hfc2pms(xyz2hfc(...)) versus xyz2pms.

## Physical / mathematical content

The test uses a random symmetric susceptibility tensor and random coordinates. Its second comparison uses the 15N isotope and checks the hyperfine-to-paramagnetic-shift tensor conversion.

## Numerical / algorithmic content

The scalar comparison passes when the absolute difference is below 1e-6; the tensor comparison uses a Frobenius-norm difference below 1e-6.

## Implementation structure

- Generate a symmetric random susceptibility tensor and compare ppcs with the isotropic trace of xyz2pms.
- Generate random coordinates, evaluate xyz2hfc followed by hfc2pms for 15N, and compare the resulting tensor with xyz2pms.
