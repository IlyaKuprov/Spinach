% A loose implementation of single-nucleon Hamiltonian eigenfunction
% calculation in the three-dimensional Woods-Saxon potential. Units,
% when not SI, are femtometres and MeV.
%
% ilya.kuprov@weizmann.ac.il

function woods_saxon(mass_number,level_number)

% Fundamental constants
r0=1.25; V0=50; a=0.5;
hbar=1.054571817e-34; % SI
pmass=1.6726219e-27;  % SI
eV=1.602176634e-19;   % SI

% Defaults
if nargin==0
    mass_number=20; level_number=5;
end

% Nuclear radius
r_nuc=r0*mass_number^(1/3);

% Simulation box dimensions
box_extents=3*[-r_nuc,r_nuc,-r_nuc,r_nuc,-r_nuc,r_nuc];
box_sizes=6*[r_nuc r_nuc r_nuc];
box_npts=[50 50 50];

% Laplacian part
premult=1e-6*1e30*hbar^2/(2*pmass*eV);
L=premult*fdlap(box_npts,box_sizes,5);

% Potential part
X=linspace(box_extents(1),box_extents(2),box_npts(1));
Y=linspace(box_extents(3),box_extents(4),box_npts(2));
Z=linspace(box_extents(5),box_extents(6),box_npts(3));
[X,Y,Z]=ndgrid(X,Y,Z); R=sqrt(X.^2+Y.^2+Z.^2);
V=-V0./(1+exp((R-r_nuc)/a));

% Plot the potential
kfigure(); volplot(V,box_extents);
ktitle(['Woods-Saxon potential, M=' num2str(mass_number)]);
kxlabel('X, fm'); kylabel('Y, fm'); kzlabel('Z, fm');

% Assemble the Hamiltonian
H=spdiags(V(:),0,prod(box_npts),prod(box_npts))-L;

% Get the state
[psi,E]=eigs(H,level_number,'smallestreal');
E=diag(E); disp('Energies, MeV:'); disp(E);

% Plot the state
psi=reshape(real(psi(:,level_number)),box_npts);
kfigure(); volplot(psi,box_extents);
ktitle(['Eig ' num2str(level_number) ' ,real part']); 
kxlabel('X, fm'); kylabel('Y, fm'); kzlabel('Z, fm');

end

