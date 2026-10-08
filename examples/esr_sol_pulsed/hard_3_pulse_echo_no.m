% DEER spin echo for a pair of nitroxide radicals at X-band. Two nit-
% roxide radicals are positioned at a distance of 25 Angstroms.
%
% The calculation is done by brute-force time propagation and numerical
% powder averaging in Liouville space. Nitroxide g-tensor data comes
% from http://dx.doi.org/10.1063/1.1697233
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function hard_3_pulse_echo_no()

% Spin system properties
sys.magnet=0.33;
sys.isotopes={'E','E'};
inter.zeeman.eigs{1}=[2.0089 2.0061 2.0027];
inter.zeeman.eigs{2}=[2.0089 2.0061 2.0027];
inter.zeeman.euler{1}=[1.0 2.0 3.0];
inter.zeeman.euler{2}=[3.0 1.0 2.0];
inter.coordinates={[ 0.00 0.00 0.00]
                   [25.00 0.00 0.00]};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.ex_prob=(operator(spin_system,{'L+'},{1})+...
                    operator(spin_system,{'L-'},{1}))/2;
parameters.ex_pump=(operator(spin_system,{'L+'},{2})+...
                    operator(spin_system,{'L-'},{2}))/2;
parameters.spins={'E'};
parameters.ta=1.0e-6;
parameters.tb=0.5e-6;
parameters.tc=0.5e-6;
parameters.nsteps=256;
parameters.grid='rep_2ang_3200pts_sph';

% Pulse sequence
echo=powder(spin_system,@deer_3p_hard_echo,parameters,'deer');

% Build the time axis
time_axis=linspace(-parameters.tc/2,...
                   +parameters.tc/2,parameters.nsteps+1);

% Plotting
kfigure(); plot(1e6*time_axis,imag(echo)); kgrid;
kxlabel('time, microseconds'); axis tight;

end

