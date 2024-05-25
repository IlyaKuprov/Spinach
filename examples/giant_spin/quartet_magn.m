% Sample magnetisation during a finite-speed magnetic field 
% sweep for a spin-3/2 particle with a zero-field splitting.
%
% Calculation time: seconds
%
% e.suturina@soton.ac.uk
% i.kuprov@soton.ac.uk

function quartet_magn()

% This must be set to 1 Tesla
sys.magnet=1.0;

% Particle
sys.isotopes={'E4'};

% Zeeman tensor
inter.zeeman.matrix={[2 0 0; 0 2 0; 0 0 2]};

% Zero-field splitting
D=icm2hz(-0.5); E=0.3*D;
inter.coupling.matrix{1,1}=zfs2mat(D,E,0,0,0);

% Formalism and basis set
bas.approximation='none';
bas.formalism='zeeman-hilb';

% Temperature
inter.temperature=1.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.fields=[0 1];
parameters.npoints=1000;
parameters.sweep_time=1e-9;      % seconds
parameters.orientation=[0 0 0];
parameters.nstates=4;

% Run the field scan
[fields,z_magn]=fieldscan_magn(spin_system,parameters);

% Plot the results
figure(); plot(fields,z_magn); kgrid;
kxlabel('Magnetic field, Tesla');
kylabel('Sample magnetisation');

end

