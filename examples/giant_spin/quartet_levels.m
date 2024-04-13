% Energy levels magnetic field scan for a spin-3/2 particle
% with a zero-field splitting.
%
% Calculation time: seconds
%
% e.suturina@soton.ac.uk
% i.kuprov@soton.ac.uk

function quartet_levels()

% This must be set to 1 Tesla
sys.magnet=1.0;

% Particle
sys.isotopes={'E4'};

% Zeeman tensor
inter.zeeman.matrix={[2 0 0; 0 2 0; 0 0 2]};

% Zero-field splitting
D=icm2hz(-0.5); E=0.3*D;
inter.coupling.matrix{1,1}=zfs2mat(D,E);

% Formalism and basis set
bas.approximation='none';
bas.formalism='zeeman-hilb';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.fields=[0 1];
parameters.npoints=100;
parameters.orientation=[0 0 0];
parameters.nstates=4;

% Run the field scan
fieldscan_enlev(spin_system,parameters)

end

