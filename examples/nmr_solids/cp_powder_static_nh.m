% 1H-15N cross-polarisation experiment in the doubly rotating
% frame. Static powder simulation.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function cp_powder_static_nh()

% System specification
sys.magnet=9.394;
sys.isotopes={'15N','1H'};
          
% Interactions
inter.zeeman.scalar={0.00 0.00};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.05]};
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
parameters.irr_opers={operator(spin_system,'Ly','1H') ...
                      operator(spin_system,'Lx','15N')};
parameters.exc_opers={operator(spin_system,'Lx','1H') ...
                      operator(spin_system,'Ly','15N')};
parameters.coil=state(spin_system,'Lx','15N');
parameters.grid='rep_2ang_6400pts_sph';
parameters.time_steps=1e-5*ones(1,100);
parameters.needs={'aniso_eq'};
parameters.spins={'15N'};

% Simulation
fid=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Plotting
time_axis=[0 cumsum(parameters.time_steps)];
figure(); plot(time_axis,real(fid)); kgrid;
kylabel('$S_{\rm{X}}$ expectation value on $^{15}N$');  
kxlabel('time, seconds'); xlim tight;

end

