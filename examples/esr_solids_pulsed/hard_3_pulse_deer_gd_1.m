% Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to
% reproduce Figure 2b from the paper by Otting and co-authors:
%
%                 http://dx.doi.org/10.1021/ja204415w
%
% The calculation is done by brute-force time propagation and grid pow-
% der averaging. Central transitions are used on both gadolinium ions.
%
% Note: gadolinium spin echo is very sharp and difficult to catch in
%       simulations because they do not include zero-field splitting
%       distributions found in experimental systems.
%
% Note: flip-flop terms in the inter-electron dipolar interaction are
%       switched off ('deer-zz') to mimic the effects of slightly dif-
%       ferent pulse frequencies in the experiment.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
% nurit.manukovsky@weizmann.ac.il
% daniella.goldfarb@weizmann.ac.il

function hard_3_pulse_deer_gd_1()

% Spin system properties
sys.magnet=3.5;
sys.isotopes={'E8','E8'};
inter.zeeman.scalar{1}=2.002319;
inter.zeeman.scalar{2}=2.002319;
inter.coordinates={[ 0.00 0.00 0.00]
                   [60.50 0.00 0.00]};
inter.coupling.matrix{1,1}=[1e8  0   0
                             0  1e8  0
                             0   0 -2e8];
inter.coupling.matrix{2,2}=[1e8  0   0
                             0  1e8  0
                             0   0 -2e8];

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.rho0=state(spin_system,'Lz','E8');
parameters.ex_prob=operator(spin_system,'CTx',1);
parameters.ex_pump=operator(spin_system,'CTx',2);
parameters.coil_prob=state(spin_system,{'L+'},{1});
parameters.coil_pump=state(spin_system,{'L+'},{2});
parameters.spectrum_sweep=1e10;
parameters.spectrum_nsteps=1024;
parameters.ex_hard=operator(spin_system,'Lx','electrons');
parameters.stepsize=1e-7;
parameters.nsteps=80;
parameters.spins={'E8'};
parameters.grid='rep_2ang_1600pts_sph';
parameters.output='detailed';

% Pulse sequence
deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer-zz');

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

