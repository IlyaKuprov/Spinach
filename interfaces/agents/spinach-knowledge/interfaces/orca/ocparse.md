# interfaces/orca/ocparse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/orca/ocparse.m`
- Signature: `[density,ext,dx,dy,dz]=ocparse(filename,pad_factor)`
- Total lines: 95

## Purpose

ORCA cube file parser. Extracts the normalised probability density and the associated metric information from ORCA spin density in "3D simple format" (see ORCA manual). Syntax: [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)

## Physical / mathematical content

- ORCA interfaces. They recover quantum-chemistry tensors and metadata and convert them to Spinach conventions.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(filename,pad_factor)`.
- Lines 36-37: Inform the user; implemented by `disp('Parsing ORCA cube ')`.
- Lines 39-40: Parse the file; implemented by `A=importdata(filename,' ',4)`.
- Lines 42-43: Get the number of points along [x y z]; implemented by `npts=str2num(A.textdata{2})`.
- Lines 45-46: Make the probability density cube; implemented by `density=reshape(abs(A.data),npts(2),npts(3),npts(1))`.
- Lines 49-50: Get the corner coordinates; implemented by `corner_xyz=str2num(A.textdata{3})`.
- Lines 52-53: Get the grid spacing vector; implemented by `dxdydz=str2num(A.textdata{4})`.
- Lines 56-57: Integrate and normalise probability density; implemented by `total_prob=trapz(trapz(trapz(density)))*dx*dy*dz`.
- Lines 60-63: Compute grid extents; implemented by `ext=[corner_xyz(1) (corner_xyz(1)+(npts(1)-1)*dx) corner_xyz(2) (corner_xyz(2)+(npts(2)-1)*dy) corner_xyz(3) (corner_xyz(3)+(npts(3)-1)*dz)]`.
- Lines 65-66: Pad the density with zeros; implemented by `density=padarray(density,pad_factor*size(density),0,'both')`.

### Key state/data transformations

- Lines 40: computes `A` using `A=importdata(filename,' ',4)`.
- Lines 43: computes `npts` using `npts=str2num(A.textdata{2})`.
- Lines 46: computes `density` using `density=reshape(abs(A.data),npts(2),npts(3),npts(1))`.
- Lines 50: computes `corner_xyz` using `corner_xyz=str2num(A.textdata{3})`.
- Lines 53: computes `dxdydz` using `dxdydz=str2num(A.textdata{4})`.
- Lines 54: computes `dx` using `dx=dxdydz(1); dy=dxdydz(2); dz=dxdydz(3)`.
- Lines 57: computes `total_prob` using `total_prob=trapz(trapz(trapz(density)))*dx*dy*dz`.
- Lines 61-63: computes `ext` using `ext=[corner_xyz(1) (corner_xyz(1)+(npts(1)-1)*dx) corner_xyz(2) (corner_xyz(2)+(npts(2)-1)*dy) corner_xyz(3) (corner_xyz(3)+(npts(3)-1)*dz)]`.

### Local helper functions

- Line 74: `grumble()` — `function grumble(filename,pad_factor)`.
  - Representative operation: `if (~ischar(filename))||isempty(filename)`.
  - Representative operation: `error('filename must be a non-empty character string.')`.

## Parameters / inputs

- filename -character string specifying the file to load
- pad_factor -padding factor specifying how many multiples
- the array dimension in zeros to add on each
- side of the cube

## Outputs

- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]
- dx,dy,dz -grid steps in the three directions, Angstrom

## Implementation structure

- ORCA cube file parser. Extracts the normalised probability density and
- the associated metric information from ORCA spin density in "3D simple
- format" (see ORCA manual). Syntax:
- [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)
- filename -character string specifying the file to load
- pad_factor -padding factor specifying how many multiples
- the array dimension in zeros to add on each
- side of the cube
- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `importdata()`, `str2num()`, `npts()`, `dxdydz()`, `trapz()`, `corner_xyz()`, `padarray()`, `ext()`, `ischar()`, `exist()`, `isscalar()`.
