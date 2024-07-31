% Gadolinium(III) DEER echo experiment. The calculation is done by
% brute-force time propagation and powder averaging. Outermost ZFS
% transition is excited by the probe pulse and the central transi-
% tion is excited by the pump pulse. Pulses are assumed to be hard.
%
% Note: gadolinium spin echo is very sharp and difficult to catch in
%       simulations because they do not include zero-field splitting
%       distributions found in experimental systems.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% nurit.manukovsky@weizmann.ac.il

function hard_3_pulse_echo_gd()

% Spin system properties
sys.magnet=3.5; D=0.56e9;
sys.isotopes={'E8','E8'};
inter.zeeman.scalar={2.002319 2.002319};
inter.coordinates={[0 0 0]; 29.5*[0 1 0]};
inter.coupling.eigs{1,1}=[D D -2*D]/3;
inter.coupling.eigs{2,2}=[D D -2*D]/3;
inter.coupling.euler{1,1}=[0 0 0];
inter.coupling.euler{2,2}=[0 pi/2 0];

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';
               
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Probe pulse operator
sigma=pauli(2); sigma.p=[zeros(6,8); [zeros(2,6) sigma.p]];
Ep_prob=kron(sigma.p,speye(size(sigma.p)));

% Pump pulse operator
sigma=pauli(2); sigma.p=[zeros(3,8); [zeros(2,3) sigma.p zeros(2,3)]; zeros(3,8)];
Ep_pump=kron(speye(size(sigma.p)),sigma.p);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E8');
parameters.ex_prob=(Ep_prob+Ep_prob')/2; 
parameters.ex_pump=(Ep_pump+Ep_pump')/2;
parameters.coil=state(spin_system,{'L+'},{1});
parameters.spins={'E8'};
parameters.ta=2e-6;
parameters.tb=1e-6;
parameters.tc=0.5e-7;
parameters.nsteps=500;
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

