% E-15N cross-polarization experiment in the doubly rotating
% frame. Static powder simulation.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function crosspol_powder_static_1()

% System specification
sys.magnet=9.394;
sys.isotopes={'15N','E'};
          
% Interactions
inter.zeeman.scalar={0.00 2.0023193043622};
inter.coordinates={[0.00 0.00 0.00 ]
                   [0.00 0.00 10.05]};
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.irr_powers=[5e4*ones(1,100);
                       5e4*ones(1,100)];
parameters.irr_opers={operator(spin_system,'Ly','E'), ...
                      operator(spin_system,'Lx','15N')};
parameters.exc_opers={operator(spin_system,'Lx','E'), ...
                      operator(spin_system,'Ly','15N')};
parameters.needs={'iso_eq'}; % Good enough here
parameters.coil=state(spin_system,'Lx','15N');
parameters.grid='rep_2ang_6400pts_sph';
parameters.time_steps=1e-5*ones(1,100);
parameters.spins={'15N'};

% Simulation
fid=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Time axis generation
t_axis=[0 cumsum(parameters.time_steps)];

% Plotting
figure(); plot(t_axis,real(fid)); kgrid;
kylabel('$S_\textrm{x}$ expectation value on $^{15}$N'); 
kxlabel('Contact pulse duration, seconds'); xlim tight;

end

