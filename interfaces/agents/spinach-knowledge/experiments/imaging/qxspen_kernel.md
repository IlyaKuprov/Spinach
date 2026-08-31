# experiments/imaging/qxspen_kernel.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/qxspen_kernel.m`
- Signature: `[K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)`
- Total lines: 104

## Purpose

Distortion kernel of the QxSPEN experiment and its derivatives. Syntax: [K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(FOVy,NSR,Nyacq,alp,bet,gam,del)`.
- Lines 46-47: Location grid in the [...] dimension; implemented by `yacq=linspace(-FOVy/2,FOVy/2,Nyacq)'`.
- Lines 49-50: Location grid and step in the [...] dimension; implemented by `ySR=linspace(-FOVy/2,FOVy/2,NSR); deltaySR=ySR(2)-ySR(1)`.
- Lines 52-53: Preallocate the answers; implemented by `K=zeros([Nyacq NSR],'like',1i)`.
- Lines 55-56: Loop over [...]; implemented by `for m=1:NSR`.
- Lines 58-60: Location grid in the integration dimension; implemented by `yint=linspace(ySR(m)-deltaySR/2, ySR(m)+deltaySR/2,100)`.
- Lines 62-63: Compute kernel element; implemented by `K(:,m)=trapz(sinc(alp*(yint-yacq)).*`.
- Lines 71-72: Compute the derivatives; implemented by `dK_dgam=1i*yacq.*K; dK_ddel=1i*K`.
- Lines 74-75: Normalise the outputs; implemented by `nfactor=norm(K,2); K=K/nfactor`.

### Control flow inferred from the code

- Line 56: `for` loop over `m=1:NSR`.

### Key state/data transformations

- Lines 47: computes `yacq` using `yacq=linspace(-FOVy/2,FOVy/2,Nyacq)'`.
- Lines 50: computes `ySR` using `ySR=linspace(-FOVy/2,FOVy/2,NSR); deltaySR=ySR(2)-ySR(1)`.
- Lines 53: computes `K` using `K=zeros([Nyacq NSR],'like',1i)`.
- Lines 59-60: computes `yint` using `yint=linspace(ySR(m)-deltaySR/2, ySR(m)+deltaySR/2,100)`.
- Lines 63: computes `K(:,m)` using `K(:,m)=trapz(sinc(alp*(yint-yacq)).*`.
- Lines 72: computes `dK_dgam` using `dK_dgam=1i*yacq.*K; dK_ddel=1i*K`.
- Lines 75: computes `nfactor` using `nfactor=norm(K,2); K=K/nfactor`.
- Lines 77: computes `dK_ddel` using `dK_ddel=dK_ddel/nfactor`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(FOVy,NSR,Nyacq,alp,bet,gam,del)`.
  - Representative operation: `if (~isnumeric(FOVy))||(~isreal(FOVy))||(~isscalar(FOVy))|| (~isnumeric(NSR))||(~isreal(NSR))||(~isscalar(NSR))|| (~isnumeric(Nyacq))||(~isreal(Nyacq))||(~isscalar(Nyacq…`.
  - Representative operation: `(~isnumeric(NSR))||(~isreal(NSR))||(~isscalar(NSR))|| (~isnumeric(Nyacq))||(~isreal(Nyacq))||(~isscalar(Nyacq))|| (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))|| (…`.

## Parameters / inputs

- FOVy -field of view along Y, mm
- NSR -number of points in the reconstruction
- by regularisation
- Nyacq -number of points that is acquired
- by the instrument
- alp -an uncertain magic number that IK does
- not understand, ask Ke Dai
- bet -an uncertain magic number that IK does
- not understand, ask Ke Dai
- gam -slope of the linear phase, rad/mm
- del -constant phase, radians

## Outputs

- K -QxSPEN kernel matrix
- dK_dgam -derivative of K with respect to gam
- dK_ddel -derivative of K with respect to del
- Note: a dangerous numerical integration stage with a fixed point count
- is used -must be replaced with an analytical expression!

## Implementation structure

- Distortion kernel of the QxSPEN experiment and its derivatives. Syntax:
- [K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)
- FOVy -field of view along Y, mm
- NSR -number of points in the reconstruction
- by regularisation
- Nyacq -number of points that is acquired
- by the instrument
- alp -an uncertain magic number that IK does
- not understand, ask Ke Dai
- bet -an uncertain magic number that IK does
- gam -slope of the linear phase, rad/mm
- del -constant phase, radians

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ySR()`, `trapz()`, `sinc()`, `isscalar()`.
