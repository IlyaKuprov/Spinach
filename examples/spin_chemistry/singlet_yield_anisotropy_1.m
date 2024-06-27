% Singlet yield anisotropy calculation for a radical pair
% using exponential recombination kinetics model.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk

function singlet_yield_anisotropy_1()

% Unit magnet (field sweep)
sys.magnet=1;

% System specification
sys.isotopes={'E','E','1H'};
inter.zeeman.scalar={2.0023 2.0023 0};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,3}=gauss2mhz([5e6   0    0; 
                                      0   4e6    0;
                                      0     0 10e6]);
% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.fields=50e-6;
parameters.rates=2e6;
parameters.electrons=[1 2];
parameters.grid='leb_2ang_rank_71';
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

