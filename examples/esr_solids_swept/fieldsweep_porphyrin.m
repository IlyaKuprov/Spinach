% Field swept EPR spectrum of copper porphyrin complex, computed
% by finding resonance fields and transition moments.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk

function fieldsweep_porphyrin()

% Magnet field
sys.magnet=1;

% Isotopes
sys.isotopes={'14N','14N','14N','14N','E','63Cu'};

% Array preallocation
inter.zeeman.eigs=cell(6,1);
inter.zeeman.euler=cell(6,1);
inter.coupling.eigs=cell(6,6);
inter.coupling.euler=cell(6,6);
inter.coupling.scalar=cell(6,6);

% Zeeman interactions
inter.zeeman.eigs{5,1}=[2.0509  2.0509  2.1801];
inter.zeeman.euler{5,1}=[0 0 0];

% Hyperfine interactions
inter.coupling.eigs{5,6}=[-70.9257  -70.9257  -575.0219]*1e6;
inter.coupling.euler{5,6}=[0 0 0];
inter.coupling.scalar{1,5}=46.0345*1e6;
inter.coupling.scalar{2,5}=46.0345*1e6;
inter.coupling.scalar{3,5}=46.0345*1e6;
inter.coupling.scalar{4,5}=46.0345*1e6;
                        
% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Symmetry
bas.sym_group={'S4'};
bas.sym_spins={[1 2 3 4]};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E'};
parameters.grid=6;
parameters.mw_freq=9.39e9;
parameters.fwhm=5e-4;
parameters.int_tol=1;
parameters.tm_tol=0.1;
parameters.window=[0.27 0.35];
parameters.npoints=512;
parameters.rspt_order=Inf;

% Run the simulation
parameters.rho0=state(spin_system,'Lz','E');
[b_axis,spec]=fieldsweep(spin_system,parameters);

% Plotting
figure(); plot(b_axis',spec');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
axis tight; kgrid;

end

