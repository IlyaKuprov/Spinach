% Gadolinium(III) DEER experiment. The calculation is done by brute-
% force time propagation and powder averaging. Outermost ZFS transi-
% tion is excited by the probe pulse and the central transition is
% excited by the pump pulse. The pulses are assumed to be ideal.
%
% Note: gadolinium spin echo is very sharp and difficult to catch in
%       simulations because they do not include zero-field splitting
%       distributions found in experimental systems.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
% nurit.manukovsky@weizmann.ac.il
% daniella.goldfarb@weizmann.ac.il

function hard_3_pulse_deer_gd_2()

% Spin system properties
sys.magnet=3.5;
sys.isotopes={'E8','E8'};
inter.zeeman.scalar={2.002319 2.002319};
inter.coupling.eigs{1,1}=[1e9 1e9 -2e9];
inter.coupling.euler{1,1}=[0 pi/9 0];
inter.coupling.eigs{2,2}=[1e9 1e9 -2e9];
inter.coupling.euler{2,2}=[0 4*pi/9 0];
inter.coordinates={[0 0 0]; 30*[sind(20) 0 cosd(20)]};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';
               
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Probe pulse operator - bottom transition
sigma=pauli(2); sigma.p=[zeros(6,8); [zeros(2,6) sigma.p]];
Ep_prob=kron(sigma.p,speye(size(sigma.p)));

% Pump pulse operator - central transition
sigma=pauli(2); sigma.p=[zeros(3,8); [zeros(2,3) sigma.p zeros(2,3)]; zeros(3,8)];
Ep_pump=kron(speye(size(sigma.p)),sigma.p);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E8');
parameters.ex_prob=(Ep_prob+Ep_prob')/2; 
parameters.ex_pump=(Ep_pump+Ep_pump')/2;
parameters.coil_prob=state(spin_system,{'L+'},{1});
parameters.coil_pump=state(spin_system,{'L+'},{2});
parameters.spectrum_sweep=5e10;
parameters.spectrum_nsteps=1024;
parameters.ex_hard=operator(spin_system,'Lx','electrons');
parameters.spins={'E8'};
parameters.stepsize=2e-8;
parameters.nsteps=100;
parameters.grid='rep_2ang_1600pts_sph';
parameters.output='detailed';

% Pulse sequence
deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer');

% Apodisation
deer.hard_pulse_fid=apodisation(spin_system,deer.hard_pulse_fid,{{'exp',6}});
deer.pump_pulse_fid=apodisation(spin_system,deer.pump_pulse_fid,{{'exp',6}});
deer.prob_pulse_fid=apodisation(spin_system,deer.prob_pulse_fid,{{'exp',6}});

% Fourier transforms
hard_pulse_spec=imag(fftshift(fft(deer.hard_pulse_fid,4*parameters.spectrum_nsteps)));
pump_pulse_spec=imag(fftshift(fft(deer.prob_pulse_fid,4*parameters.spectrum_nsteps)));
prob_pulse_spec=imag(fftshift(fft(deer.pump_pulse_fid,4*parameters.spectrum_nsteps)));
freq_axis_hz=ft_axis(0,parameters.spectrum_sweep,4*parameters.spectrum_nsteps);

% Plotting
figure(); scale_figure([1.0 2.0]);
time_axis=(0:parameters.nsteps)*parameters.stepsize;
subplot(4,1,1); plot(freq_axis_hz,hard_pulse_spec);
kxlabel('Offset frequency, Hz'); xlim tight; kgrid; 
title('frequency swept spectrum');
subplot(4,1,2); plot(freq_axis_hz,pump_pulse_spec); 
kxlabel('Offset frequency, Hz'); xlim tight; kgrid; 
title('excitation profile, probe spin');
subplot(4,1,3); plot(freq_axis_hz,prob_pulse_spec); 
kxlabel('Offset frequency, Hz'); xlim tight; kgrid; 
title('excitation profile, pump spin');
subplot(4,1,4); plot(time_axis,-imag(deer.deer_trace));
kxlabel('time, seconds'); xlim tight; kgrid; 
ktitle('DEER trace');

end

