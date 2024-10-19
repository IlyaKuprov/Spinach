% A comparison between a point model fit and a multipole model 
% fit described in
%
%            http://dx.doi.org/10.1039/c6cp05437d
%
% in a situation where the probability density of the paramag-
% netic centre is certainly not a point (four randomly placed
% Gaussians are used).
%
% Calculation time: minutes
%
% e.suturina@soton.ac.uk

function point_vs_distr()

% Size of the cube containing the Gaussians
bounding_cube_size=3; % Angstrom

% Randomly positioned Gaussian centroids
c=bounding_cube_size*(rand(4,3))-bounding_cube_size/2;

% Gaussian standard deviations
sigma=0.5; % Angstrom

% number of points in each dimension of the 3D grid
npoints=64;

% Grid extents
blind_radius=5; layer_thickness=10; % Angstrom
a=blind_radius+layer_thickness;
ext=[-a a -a a -a a];

% Get a 3D grid
[X,Y,Z]=ndgrid(linspace(-a,a,npoints),...
               linspace(-a,a,npoints),...
               linspace(-a,a,npoints));

% Get the normalized density
rho=1/(sqrt(2*pi)*sigma)^3*(exp(-((X-c(1,1)).^2+(Y-c(1,2)).^2+(Z-c(1,3)).^2)/(2*sigma^2))+...
                            exp(-((X-c(2,1)).^2+(Y-c(2,2)).^2+(Z-c(2,3)).^2)/(2*sigma^2))+...
                            exp(-((X-c(3,1)).^2+(Y-c(3,2)).^2+(Z-c(3,3)).^2)/(2*sigma^2))+...
                            exp(-((X-c(4,1)).^2+(Y-c(4,2)).^2+(Z-c(4,3)).^2)/(2*sigma^2)))/4;
                        
% Number of randomly placed nuclei
n_nuclei=100;

% Random positions for the nuclei within a spherical layer
theta=pi*rand(n_nuclei,1);
phi=2*pi*rand(n_nuclei,1);
blind_radius=5; layer_thickness=10;
r=layer_thickness*rand(n_nuclei,1)+blind_radius;
x=r.*sin(theta).*cos(phi);
y=r.*sin(theta).*sin(phi);
z=r.*cos(theta);

% Plot density and nuclei
volplot(rho,ext); hold on; 
plot3(x,y,z,'b.'); drawnow();

% A random susteptibility tensor (cubic Angstroms)
ax=-0.45; rh=-0.05; R=euler2dcm(pi/3,pi/4,pi/5);
chi=R*diag([(-ax/3+rh) (-ax/3-rh) (2*ax/3)])*R';

% Pad density with zeros (two volumes each side) to avoid PBC effects
pad_size=2; pad_density=padarray(rho,pad_size*size(rho),0,'both');

% Update grid extents
ext=(2*pad_size+1)*ext;

% Compute PCS exactly using FFT method for PDE equation
[expt_pcs,~]=kpcs(pad_density,chi,ext,[x y z],'fft');

% Fit the results with the point model
[mxyz_p,chi_p,theo_pcs_p,s_mxyz_p,s_chi_p]=ippcs([x y z],[0 0 0],expt_pcs);

% Define the ranks of multipole moments to be considered
ranks=[0 1 2];

% Fit the results with the multipole model up to second rank
[mxyz_m,chi_m,Ilm,theo_pcs_m,s_mxyz_m,s_chi_m,s_Ilm]=ilpcs([x y z],expt_pcs,ranks,[0 0 0]);

% Print the true values
disp(' '); disp('TRUE VALUES'); disp('-----------');
disp('Susceptibility tensor (cubic Angstrom):'); disp(chi);

% Print the point model results
disp(' '); disp('POINT APPROXIMATION');
disp('-------------------');
disp('Paramagnetic centre position (Angstrom):'); disp(mxyz_p);
disp('Standard deviation (Angstrom):');           disp(s_mxyz_p);
disp(' ');
disp('Susceptibility tensor (cubic Angstrom):');  disp(chi_p);
disp('Standard deviation (cubic Angstrom):   ');  disp(s_chi_p);

% Print the multipole model results
disp(' '); disp('MULTIPOLE APPROXIMATION');
disp('-----------------------');
disp('Origin (Angstrom):'); disp(mxyz_m);
disp('Standard deviation (Angstrom):');           disp(s_mxyz_m);
disp(' ');
disp('Susceptibility tensor (cubic Angstrom):');  disp(chi_m);
disp('Standard deviation (cubic Angstrom):   ');  disp(s_chi_m);

% % Plot the diagnostics
figure();   plot(sqrt(x.^2+y.^2+z.^2),theo_pcs_p-expt_pcs,'bo');
hold('on'); plot(sqrt(x.^2+y.^2+z.^2),theo_pcs_m-expt_pcs,'ro'); kgrid;
kxlabel('Distance from the origin, Angstrom'); kylabel('PCS fitting error (ppm)')
klegend('point model','delocalised model','Location','northeast'); hold('off');

% Compare fitted multipole moments with the original ones
Ilm0=points2mult([X(:) Y(:) Z(:)],mxyz_m,rho(:),[0 1 2],'grid');
disp('Moments of the true density:'); disp(cell2mat(Ilm0));
disp('Moments recovered from PCS fitting:'); disp(cell2mat(Ilm));
disp('Standard deviations from PCS fitting:'); disp(cell2mat(s_Ilm));

end

