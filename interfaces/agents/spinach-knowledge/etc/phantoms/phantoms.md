# etc/phantoms/phantoms.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/phantoms/phantoms.m`
- Signature: `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)`
- Total lines: 169

## Purpose

MRI phantom library. Syntax: [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(ph_name)`.
- Lines 31-32: Get own location; implemented by `P=mfilename('fullpath'); P=P(1:(end-9))`.
- Lines 38-39: Load the breast phantom; implemented by `load([P filesep 'mri_lab_boob.mat'],'boob','nvoxels')`.
- Lines 41-42: Preallocate arrays; implemented by `R1Ph=zeros(nvoxels); R2Ph=zeros(nvoxels); PDPh=zeros(nvoxels)`.
- Lines 44-45: Set number of points; implemented by `npts=size(R1Ph)`.
- Lines 47-48: Assume 1 mm resolution; implemented by `dims=npts.*[1e-3 1e-3 1e-3]`.
- Lines 50-51: Assign water proton density (1.0 for pure water); implemented by `PDPh(boob==-1.0)=0.00`.
- Lines 62-63: Assign longitudinal water relaxation rate (3 Tesla); implemented by `R1Ph(boob==-1.0)=0.00`.
- Lines 74-75: Assign transverse water relaxation rate (3 Tesla); implemented by `R2Ph(boob==-1.0)= 0.00`.
- Lines 88-89: Load MRiLab phantom; implemented by `load([P filesep 'mri_lab_brain.mat'],'VObj')`.
- Lines 91-92: Assign the variables; implemented by `R1Ph=1./VObj.T1; R2Ph=1./VObj.T2; PDPh=VObj.Rho`.
- Lines 94-95: Clean up infinities; implemented by `R1Ph(isinf(R1Ph))=0; R2Ph(isinf(R2Ph))=0`.
- Lines 97-98: Point counts; implemented by `npts=size(PDPh)`.
- Lines 100-101: Physical dimensions; implemented by `dims=npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes]`.
- Lines 114-115: Downsample; implemented by `R1Ph=R1Ph(1:2:end,1:2:end,1:2:end)`.
- Lines 122-123: Physical dimensions; implemented by `dims=2*npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes]`.
- Lines 136-137: Downsample; implemented by `R1Ph=R1Ph(1:4:end,1:4:end,1:4:end)`.
- Lines 144-145: Physical dimensions; implemented by `dims=4*npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes]`.

### Control flow inferred from the code

- Line 34: dispatches on `ph_name`; cases `'boob'`, `'brain-highres'`, `'brain-medres'`, `'brain-lowres'`.

### Key state/data transformations

- Lines 32: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-9))`.
- Lines 42: computes `R1Ph` using `R1Ph=zeros(nvoxels); R2Ph=zeros(nvoxels); PDPh=zeros(nvoxels)`.
- Lines 45: computes `npts` using `npts=size(R1Ph)`.
- Lines 48: computes `dims` using `dims=npts.*[1e-3 1e-3 1e-3]`.
- Lines 95: computes `R1Ph(isinf(R1Ph))` using `R1Ph(isinf(R1Ph))=0; R2Ph(isinf(R2Ph))=0`.
- Lines 116: computes `R2Ph` using `R2Ph=R2Ph(1:2:end,1:2:end,1:2:end)`.
- Lines 117: computes `PDPh` using `PDPh=PDPh(1:2:end,1:2:end,1:2:end)`.

### Local helper functions

- Line 157: `grumble()` — `function grumble(ph_name)`. Thank you for saving me from the Library. I can see a lot of girls without boyfriends there, and all they
  - Representative operation: `if ~ischar(ph_name)`.
  - Representative operation: `error('ph_name must be a character string.')`.

## Parameters / inputs

- ph_name -character string giving the name of the
- phantom (see the function text)

## Outputs

- R1Ph -a cube of R1 values
- R2Ph -a cube of R2 values
- PDPh -a cube of PD values
- dims -row vector of three cube dimensions, m
- npts -row vector of three cube dimensions, points

## Implementation structure

- MRI phantom library. Syntax:
- [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)
- ph_name -character string giving the name of the
- phantom (see the function text)
- R1Ph -a cube of R1 values
- R2Ph -a cube of R2 values
- PDPh -a cube of PD values
- dims -row vector of three cube dimensions, m
- npts -row vector of three cube dimensions, points
- Check consistency
- Get own location
- Load the breast phantom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `load()`, `PDPh()`, `R1Ph()`, `R2Ph()`, `isinf()`, `ischar()`.
