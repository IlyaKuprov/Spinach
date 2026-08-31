# examples/nmr_paramag/point_vs_distr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/point_vs_distr.m`
- Signature: `point_vs_distr()`
- Total lines: 117

## Purpose

A comparison between a point model fit and a multipole model fit described in in a situation where the probability density of the paramag- netic centre is certainly not a point (four randomly placed Gaussians are used). Calculation time: minutes

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Size of the cube containing the Gaussians; implemented by `bounding_cube_size=3`.
- Lines 19-20: Randomly positioned Gaussian centroids; implemented by `c=bounding_cube_size*(rand(4,3))-bounding_cube_size/2`.
- Lines 22-23: Gaussian standard deviations; implemented by `sigma=0.5`.
- Lines 25-26: number of points in each dimension of the 3D grid; implemented by `npoints=64`.
- Lines 28-29: Grid extents; implemented by `blind_radius=5; layer_thickness=10`.
- Lines 33-36: Get a 3D grid; implemented by `[X,Y,Z]=ndgrid(linspace(-a,a,npoints), linspace(-a,a,npoints), linspace(-a,a,npoints))`.
- Lines 38-42: Get the normalized density; implemented by `rho=1/(sqrt(2*pi)*sigma)^3*(exp(-((X-c(1,1)).^2+(Y-c(1,2)).^2+(Z-c(1,3)).^2)/(2*sigma^2))+ exp(-((X-c(2,1)).^2+(Y-c(2,2)).^2+(Z-c(2,3)).^2)/(2*sigma^2))+ exp(-((X-c(3,1)…`.
- Lines 44-45: Number of randomly placed nuclei; implemented by `n_nuclei=100`.
- Lines 47-48: Random positions for the nuclei within a spherical layer; implemented by `theta=pi*rand(n_nuclei,1)`.
- Lines 56-57: Plot density and nuclei; implemented by `volplot(rho,ext); hold on`.
- Lines 60-61: A random susteptibility tensor (cubic Angstroms); implemented by `ax=-0.45; rh=-0.05; R=euler2dcm(pi/3,pi/4,pi/5)`.
- Lines 64-65: Pad density with zeros (two volumes each side) to avoid PBC effects; implemented by `pad_size=2; pad_density=padarray(rho,pad_size*size(rho),0,'both')`.
- Lines 67-68: Update grid extents; implemented by `ext=(2*pad_size+1)*ext`.
- Lines 70-71: Compute PCS exactly using FFT method for PDE equation; implemented by `[expt_pcs,~]=kpcs(pad_density,chi,ext,[x y z],'fft')`.
- Lines 73-74: Fit the results with the point model; implemented by `[mxyz_p,chi_p,theo_pcs_p,s_mxyz_p,s_chi_p]=ippcs([x y z],[0 0 0],expt_pcs)`.
- Lines 76-77: Define the ranks of multipole moments to be considered; implemented by `ranks=[0 1 2]`.
- Lines 79-80: Fit the results with the multipole model up to second rank; implemented by `[mxyz_m,chi_m,Ilm,theo_pcs_m,s_mxyz_m,s_chi_m,s_Ilm]=ilpcs([x y z],expt_pcs,ranks,[0 0 0])`.
- Lines 82-83: Print the true values; implemented by `disp(' '); disp('TRUE VALUES'); disp('-----------')`.

### Key state/data transformations

- Lines 17: computes `bounding_cube_size` using `bounding_cube_size=3`.
- Lines 20: computes `c` using `c=bounding_cube_size*(rand(4,3))-bounding_cube_size/2`.
- Lines 23: computes `sigma` using `sigma=0.5`.
- Lines 26: computes `npoints` using `npoints=64`.
- Lines 29: computes `blind_radius` using `blind_radius=5; layer_thickness=10`.
- Lines 30: computes `a` using `a=blind_radius+layer_thickness`.
- Lines 31: computes `ext` using `ext=[-a a -a a -a a]`.
- Lines 34-36: computes `[X,Y,Z]` using `[X,Y,Z]=ndgrid(linspace(-a,a,npoints), linspace(-a,a,npoints), linspace(-a,a,npoints))`.
- Lines 39-42: computes `rho` using `rho=1/(sqrt(2*pi)*sigma)^3*(exp(-((X-c(1,1)).^2+(Y-c(1,2)).^2+(Z-c(1,3)).^2)/(2*sigma^2))+ exp(-((X-c(2,1)).^2+(Y-c(2,2)).^2+(Z-c(2,3)).^2)/(2*sigma^2))+ exp(-((X-c(3,1)…`.
- Lines 45: computes `n_nuclei` using `n_nuclei=100`.
- Lines 48: computes `theta` using `theta=pi*rand(n_nuclei,1)`.
- Lines 49: computes `phi` using `phi=2*pi*rand(n_nuclei,1)`.
- Lines 51: computes `r` using `r=layer_thickness*rand(n_nuclei,1)+blind_radius`.
- Lines 52: computes `x` using `x=r.*sin(theta).*cos(phi)`.
- Lines 53: computes `y` using `y=r.*sin(theta).*sin(phi)`.
- Lines 54: computes `z` using `z=r.*cos(theta)`.
- Lines 61: computes `ax` using `ax=-0.45; rh=-0.05; R=euler2dcm(pi/3,pi/4,pi/5)`.
- Lines 62: computes `chi` using `chi=R*diag([(-ax/3+rh) (-ax/3-rh) (2*ax/3)])*R'`.

## Implementation structure

- A comparison between a point model fit and a multipole model
- fit described in
- in a situation where the probability density of the paramag-
- netic centre is certainly not a point (four randomly placed
- Gaussians are used).
- Calculation time: minutes
- Size of the cube containing the Gaussians
- Randomly positioned Gaussian centroids
- Gaussian standard deviations
- number of points in each dimension of the 3D grid
- Grid extents
- Get a 3D grid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `volplot()`, `plot3()`, `euler2dcm()`, `padarray()`, `kpcs()`, `ippcs()`, `ilpcs()`, `tensor()`, `position()`, `deviation()`, `Origin()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`, `points2mult()`.
