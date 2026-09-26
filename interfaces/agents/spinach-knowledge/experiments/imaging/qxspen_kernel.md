# experiments/imaging/qxspen_kernel.m

- Signature: `[K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)`

## Purpose

Distortion kernel of the QxSPEN experiment and its derivatives. Syntax: [K,dK_dgam,dK_ddel]=qxspen_kernel(FOVy,NSR,Nyacq,alp,bet,gam,del)

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

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
