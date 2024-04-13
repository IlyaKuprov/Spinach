% Three-pulse DEER echo on a Cu(II)-NO two electron system at X-band.
%
% The calculation is done by brute-force time propagation and numerical
% powder averaging in Liouville space.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function hard_3_pulse_echo_cu()

% Spin system parameters
sys.magnet=0.33;
sys.isotopes={'E','E'};
inter.zeeman.eigs={[2.056, 2.056, 2.205];
                   [2.009, 2.006, 2.003]};
inter.zeeman.euler={[0 0 0]; [0 0 0]};
inter.coordinates={[0 0 0]; [20 0 0]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory level SSR algorithms
sys.disable={'trajlevel'};
               
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
parameters.ta=2.0e-7;
parameters.tb=1.0e-7;
parameters.tc=2.5e-8;
parameters.nsteps=256;
parameters.grid='rep_2ang_1600pts_sph';

% Pulse sequence
echo=powder(spin_system,@deer_3p_hard_echo,parameters,'deer');

% Build the time axis
time_axis=linspace(-parameters.tc/2,...
                   +parameters.tc/2,parameters.nsteps+1);

% Plotting
figure(); plot(1e6*time_axis,imag(echo)); kgrid;
xlabel('time, microseconds'); axis tight;

end

