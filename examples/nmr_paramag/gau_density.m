% Simple calculation of the PCS field of a Gaussian distribution of the
% electron probability density.
%
% g.t.p.charnock@oerc.ox.ac.uk
% ilya.kuprov@weizmann.ac.il

function gau_density()

% Set problem dimensions
dim=100; ext=10; sigma=0.5;

% Get a 3D grid
[X,Y,Z]=meshgrid(linspace(-ext,ext,dim),linspace(-ext,ext,dim),linspace(-ext,ext,dim));

% Pick a reasonable susceptibility tensor
R=euler2dcm(pi/4,pi/5,pi/6);
chi=R*[-0.1 0 0; 0 -0.2 0; 0 0 0.3]*R';

% Get electron distribution
probden=(1/sqrt((2*pi)^3*sigma^3))*exp(-(X.^2+Y.^2+Z.^2)/(2*sigma));

% Solve Kuprov equation
[~,pcs_3d]=kpcs(probden,chi,[-ext ext -ext ext -ext ext],[0 0 0],'fft');

% Plot the solution
kfigure(); volplot(pcs_3d,[-ext ext -ext ext -ext ext]); 
kgrid; ktitle('PCS field');

end

