# tests/kernel/test_transform_tensor_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_tensor_suite.m`
- Signature: `result=test_transform_tensor_suite()`
- Total lines: 149

## Purpose

Tests tensor transform helpers. Syntax: result=test_transform_tensor_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Tensor transform helpers\n')`.
- Lines 20-23: State the tensor target of the test; implemented by `result=new_test_result('kernel/transform_tensor_suite', 'Tensor transform helpers', 'tensor transforms must preserve their algebraic definitions and round-trips.')`.
- Lines 25-26: Check Haeberlen anisotropy and asymmetry at zero Euler angles; implemented by `iso=10; aniso=6; asym=0.25`.
- Lines 33-34: Check axiality and rhombicity matrix construction; implemented by `iso=4; ax=6; rh=2`.
- Lines 49-50: Check span and skew construction in the Herzfeld-Berger convention; implemented by `iso=5; span=6; skew=0.5`.
- Lines 56-57: Check zero-field splitting tensor construction; implemented by `D=9; E=2`.
- Lines 65-66: Check isotropic-antisymmetric-symmetric decomposition and reconstruction; implemented by `C=[1 2 -3;4 -5 6;7 8 9]`.
- Lines 78-79: Round-trip the Cartesian tensor through irreducible spherical components; implemented by `[rank0,rank1,rank2]=mat2sphten(C)`.
- Lines 84-85: Check spherical harmonic coefficients of an isotropic quadratic form; implemented by `[r0,r1,r2]=qform2sph(3*eye(3))`.
- Lines 93-94: Check the axial Stevens q=0 component through the rank-one transform; implemented by `Bkq=stev2sph(1,[0;1;0])`.
- Lines 98-99: Check traceless symmetric matrix parameter extraction; implemented by `T=diag([-2 -1 3])`.
- Lines 108-109: Check quadrupolar tensor construction at zero Euler angles; implemented by `Cq=1200; eta=0.2; spin_q=1`.
- Lines 115-116: Check CASTEP electric-field-gradient scaling convention; implemented by `V=diag([-1 -2 3])`.
- Lines 123-124: Check WebLab cone convention delegation for two sites; implemented by `alpha=0.11; theta=0.42; phi=0.73`.
- Lines 131-132: Check rotational averaging around the z axis; implemented by `T=diag([1 3 5])`.
- Lines 138-139: Check spin-half Hamiltonian decomposition into Zeeman terms; implemented by `S=pauli(2)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/transform_tensor_suite', 'Tensor transform helpers', 'tensor transforms must preserve their algebraic definitions and round-trips.')`.
- Lines 26: computes `iso` using `iso=10; aniso=6; asym=0.25`.
- Lines 27: computes `red_aniso` using `red_aniso=2*aniso/3`.
- Lines 28: computes `M_ref` using `M_ref=diag([iso-red_aniso*(1+asym)/2,iso-red_aniso*(1-asym)/2,iso+red_aniso])`.
- Lines 29: computes `M` using `M=anas2mat(iso,aniso,asym,0,0,0)`.
- Lines 39: computes `[iso_obs,ax_obs,rh_obs,eigs_obs]` using `[iso_obs,ax_obs,rh_obs,eigs_obs]=mat2axrh(M_ref)`.
- Lines 57: computes `D` using `D=9; E=2`.
- Lines 66: computes `C` using `C=[1 2 -3;4 -5 6;7 8 9]`.
- Lines 67: computes `[a,d,A]` using `[a,d,A]=mat2ias(C)`.
- Lines 68: computes `C_obs` using `C_obs=ias2mat(a,d,A)`.
- Lines 79: computes `[rank0,rank1,rank2]` using `[rank0,rank1,rank2]=mat2sphten(C)`.
- Lines 85: computes `[r0,r1,r2]` using `[r0,r1,r2]=qform2sph(3*eye(3))`.
- Lines 94: computes `Bkq` using `Bkq=stev2sph(1,[0;1;0])`.
- Lines 96: computes `'the q` using `'the q=0 rank-one Stevens component is already the axial spherical component')`.
- Lines 99: computes `T` using `T=diag([-2 -1 3])`.
- Lines 100: computes `[ax_obs,rh_obs,angles_obs]` using `[ax_obs,rh_obs,angles_obs]=tsm2param(T)`.
- Lines 109: computes `Cq` using `Cq=1200; eta=0.2; spin_q=1`.
- Lines 110: computes `Q_ref` using `Q_ref=diag([-240 -360 600])`.

## Outputs

- result -regression test result with explanatory messages
- The test checks interaction tensor parametrisations, spherical tensor
- round-trips, quadrupolar conversions, axial symmetrisation, and simple
- Hamiltonian decomposition.

## Implementation structure

- Tests tensor transform helpers. Syntax:
- result=test_transform_tensor_suite()
- result -regression test result with explanatory messages
- The test checks interaction tensor parametrisations, spherical tensor
- round-trips, quadrupolar conversions, axial symmetrisation, and simple
- Hamiltonian decomposition.
- Announce the test target
- State the tensor target of the test
- Check Haeberlen anisotropy and asymmetry at zero Euler angles
- Check axiality and rhombicity matrix construction
- Check span and skew construction in the Herzfeld-Berger convention
- Check zero-field splitting tensor construction

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `anas2mat()`, `test_close()`, `axrh2mat()`, `mat2axrh()`, `spsk2mat()`, `zfs2mat()`, `mat2ias()`, `ias2mat()`, `mat2sphten()`, `sphten2mat()`, `qform2sph()`, `stev2sph()`, `tsm2param()`, `euler2dcm()`, `eeqq2nqi()`.
