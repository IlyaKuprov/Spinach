# tests/kernel/test_transform_roundtrip_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_roundtrip_suite.m`
- Signature: `result=test_transform_roundtrip_suite()`
- Total lines: 85

## Purpose

Tests deterministic coordinate and tensor transforms. Syntax: result=test_transform_roundtrip_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Coordinate and tensor transform functions\n')`.
- Lines 19-22: State the transform target of the test; implemented by `result=new_test_result('kernel/transform_roundtrip_suite', 'Coordinate and tensor transform functions', 'transform helpers must preserve rotations, coordinates, and tens…`.
- Lines 24-25: Direction-cosine matrices must be orthogonal proper rotations; implemented by `R=anax2dcm([0 0 2],pi/3)`.
- Lines 31-32: Quaternion and angle-axis representations must describe the same rotation; implemented by `axis=[1 2 3]; angle=0.37*pi`.
- Lines 40-41: Euler conversion is ill-conditioned in angles, but DCM reconstruction is unique; implemented by `angles=[0.21*pi 0.37*pi 0.43*pi]`.
- Lines 47-48: Axiality/rhombicity to matrix with zero Euler angles gives the Mehring-order eigenvalues; implemented by `iso=4; ax=6; rh=2`.
- Lines 61-62: Cartesian and irreducible spherical tensor representations are algebraic inverses; implemented by `T=[1 2 3;4 5 6;7 8 10]`.
- Lines 68-69: Fractional crystallographic coordinates in an orthorhombic cell scale by cell edges; implemented by `ABC=[0 0 0; 1/2 1/3 1/4; 1 1 1]`.
- Lines 76-77: Spherical coordinates use ISO radius/inclination/azimuth convention; implemented by `[r,theta,phi]=xyz2sph([1 0 0],[0 1 0],[0 0 1])`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/transform_roundtrip_suite', 'Coordinate and tensor transform functions', 'transform helpers must preserve rotations, coordinates, and tens…`.
- Lines 25: computes `R` using `R=anax2dcm([0 0 2],pi/3)`.
- Lines 32: computes `axis` using `axis=[1 2 3]; angle=0.37*pi`.
- Lines 33: computes `q` using `q=anax2qter(axis,angle)`.
- Lines 34: computes `[axis_back,angle_back]` using `[axis_back,angle_back]=qter2anax(q)`.
- Lines 35: computes `R_from_axis` using `R_from_axis=anax2dcm(axis,angle)`.
- Lines 36: computes `R_from_quat` using `R_from_quat=anax2dcm(axis_back,angle_back)`.
- Lines 41: computes `angles` using `angles=[0.21*pi 0.37*pi 0.43*pi]`.
- Lines 43: computes `angles_back` using `angles_back=dcm2euler(R)`.
- Lines 48: computes `iso` using `iso=4; ax=6; rh=2`.
- Lines 49: computes `M` using `M=axrh2mat(iso,ax,rh,0,0,0)`.
- Lines 50: computes `M_ref` using `M_ref=diag([iso-(ax+3*rh)/6, iso-(ax-3*rh)/6, iso+ax/3])`.
- Lines 53: computes `[iso_back,ax_back,rh_back,eigvals]` using `[iso_back,ax_back,rh_back,eigvals]=mat2axrh(M_ref)`.
- Lines 62: computes `T` using `T=[1 2 3;4 5 6;7 8 10]`.
- Lines 63: computes `[r0,r1,r2]` using `[r0,r1,r2]=mat2sphten(T)`.
- Lines 64: computes `T_back` using `T_back=sphten2mat(r0,r1,r2)`.
- Lines 69: computes `ABC` using `ABC=[0 0 0; 1/2 1/3 1/4; 1 1 1]`.
- Lines 70: computes `[XYZ,va,vb,vc]` using `[XYZ,va,vb,vc]=frac2cart(2,3,4,90,90,90,ABC)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks transformation functions by using exact geometrical
- identities, algebraic inverses, and known tensor decompositions.

## Implementation structure

- Tests deterministic coordinate and tensor transforms. Syntax:
- result=test_transform_roundtrip_suite()
- result -regression test result with explanatory messages
- The test checks transformation functions by using exact geometrical
- identities, algebraic inverses, and known tensor decompositions.
- Announce the test target
- State the transform target of the test
- Direction-cosine matrices must be orthogonal proper rotations
- Quaternion and angle-axis representations must describe the same rotation
- Euler conversion is ill-conditioned in angles, but DCM reconstruction is unique
- Axiality/rhombicity to matrix with zero Euler angles gives the Mehring-order eigenvalues
- Cartesian and irreducible spherical tensor representations are algebraic inverses

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `anax2dcm()`, `test_close()`, `anax2qter()`, `qter2anax()`, `euler2dcm()`, `dcm2euler()`, `axrh2mat()`, `mat2axrh()`, `eigvals()`, `mat2sphten()`, `sphten2mat()`, `frac2cart()`, `xyz2sph()`.
