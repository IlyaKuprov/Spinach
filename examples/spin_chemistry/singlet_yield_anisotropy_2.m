% Singlet yield anisotropy calculation for a radical pair
% using exponential recombination kinetics model.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_anisotropy_2()

% Unit magnet (field sweep)
sys.magnet=1;

% Isotopes
sys.isotopes={'E','E','14N','14N','1H'};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Rotation matrices
R1=[0.4380  0.8655 -0.2432;  0.8981 -0.4097  0.1595; -0.0384  0.2883  0.9568];
R2=[0.9703 -0.2207  0.0992;  0.2383  0.9426 -0.2340; -0.0419  0.2506  0.9672];
R3=[0.9819  0.1883 -0.0203; -0.0348  0.2850  0.9579; -0.1861  0.9398 -0.2864];

% Interaction eigenvalues
A1=[-1.049  0 0; 0 -0.996 0; 0 0 13.826];
A2=[-0.305  0 0; 0 -0.222 0; 0 0 6.872];
A3=[-13.850 0 0; 0 -9.372 0; 0 0 0.143];

% Coupling tensors
inter.coupling.matrix=cell(5);
inter.coupling.matrix{1,3}=1e6*gauss2mhz(R1*A1*R1');
inter.coupling.matrix{2,4}=1e6*gauss2mhz(R2*A2*R2');
inter.coupling.matrix{1,5}=1e6*gauss2mhz(R3*A3*R3');

% Zeeman interactions
inter.zeeman.scalar={2.0023 2.0023 0 0 0};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.npoints=1;
parameters.fields=50e-6;
parameters.rates=50e6;
parameters.electrons=[1 2];
parameters.grid='leb_2ang_rank_35';
parameters.spins={'E'};
parameters.needs={'zeeman_op'};
parameters.sum_up=0;

% Simulation
[yield,grid]=powder(spin_system,@rydmr_exp,parameters,'labframe');

% Preprocessing
yield=cell2mat(yield);
yield=yield-sum(yield.*grid.weights);

% Plotting
hull=get_hull(grid.betas,grid.gammas);
x=yield.*sin(grid.betas).*cos(grid.gammas);
y=yield.*sin(grid.betas).*sin(grid.gammas);
z=yield.*cos(grid.betas); figure();
c=sqrt(x.^2+y.^2+z.^2);
trisurf(hull,x,y,z,c,'EdgeAlpha',0.25); 
kgrid; axis square; axis equal;
axis tight; box on;

end

