% 1H-15N cross-polarisation experiment in the doubly rotating
% frame using (a) tangent-ramped adiabatic CP; (b) numerically
% optimised (GRAPE method) shortcut to adiabaticity.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cp_adiabatic_vs_optimcon()

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

%% Tangent ramp CP simulation

% Common experiment parameters
parameters.time_steps=2e-6*ones(1,500);
parameters.irr_opers={operator(spin_system,'Ly','1H') ...
                      operator(spin_system,'Lx','15N')};
parameters.exc_opers={operator(spin_system,'Lx','1H')};
parameters.coil=state(spin_system,'Lx','15N');
parameters.grid='rep_2ang_200pts_sph';
parameters.needs={'aniso_eq'};
parameters.spins={'15N'};

% Simulate tangent ramped amplitude CP
ramp_up=tan(linspace(-1.4,1.4,500));
ramp_up=ramp_up-min(ramp_up); 
ramp_up=ramp_up/max(ramp_up);
ramp_down=fliplr(ramp_up);
irr_powers_a=[5e4*ramp_down; 5e4*ramp_up];
parameters.irr_powers=irr_powers_a;
fid_a=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Plotting - waveform
kfigure(); scale_figure([1.5 1.5]);
time_axis=[0 cumsum(parameters.time_steps)];
subplot(2,2,1); plot(time_axis(2:end),irr_powers_a);
xlim tight; kgrid; kxlabel('time, seconds');
kylabel('nutation frequency, Hz'); ylim([-2e4 6e4]);
klegend({'Tangent ramp $^{1}$H',...
         'Tangent ramp $^{15}$N'},'Location','South');

% Plotting - trajectory
subplot(2,2,4); plot(time_axis,real(fid_a)); hold on;
kylabel('$S_{\rm{X}}$ expectation value on $^{15}N$'); 
xlim tight; kgrid; kxlabel('time, seconds'); drawnow();

%% GRAPE optimisation with the same timing

% Drift Hamiltonians for every system in the powder
control.drifts=drifts(spin_system,@powder,parameters,'qnmr');

% Initial state: Ly on 1H
rho_init=state(spin_system,'Ly','1H');
rho_init=rho_init/norm(rho_init,2);
control.rho_init={rho_init};

% Target state: Lx on 15N
rho_targ=state(spin_system,'Lx','15N');
rho_targ=rho_targ/norm(rho_targ);
control.rho_targ={rho_targ};

% Control operators
LyH=operator(spin_system,'Ly','1H');
LxN=operator(spin_system,'Lx','15N');
control.operators={LyH,LxN};

% Other GRAPE settings
control.pwr_levels=2*pi*5e4;              % Pulse power, rad/s
control.pulse_dt=parameters.time_steps;   % Slice durations
control.penalties={'SNS'};                % Penalty
control.p_weights=100;                    % Penalty weight
control.method='lbfgs';                   % Optimisation method
control.max_iter=30;                      % Termination tolerance
control.parallel='ensemble';              % Parallelisation

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess from above
guess=[ramp_down; ramp_up];

% Run the optimisation
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Simulate GRAPE CP
parameters.irr_powers=5e4*pulse;
fid_b=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Plotting
irr_powers_b=5e4*pulse;
subplot(2,2,2); plot(time_axis(2:end),irr_powers_b);
xlim tight; kgrid; kxlabel('time, seconds');
kylabel('nutation frequency, Hz'); ylim([-2e4 6e4]);
klegend({'GRAPE (same time) $^{1}$H',...
         'GRAPE (same time) $^{15}$N'},'Location','South');
subplot(2,2,4); plot(time_axis,real(fid_b)); drawnow();

%% GRAPE optimisation with a half of the time

% Update the slice durations
parameters.time_steps=parameters.time_steps/2;
control.pulse_dt=parameters.time_steps;

% Re-run the housekeeping
spin_system=optimcon(spin_system,control);

% Re-run the optimisation
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Simulate GRAPE CP
parameters.irr_powers=5e4*pulse;
fid_c=powder(spin_system,@cp_contact_hard,parameters,'nmr');

% Plotting
irr_powers_c=5e4*pulse;
time_axis=[0 cumsum(control.pulse_dt)];
subplot(2,2,3); plot(time_axis(2:end),irr_powers_c);
xlim tight; kgrid; kxlabel('time, seconds'); 
kylabel('nutation frequency, Hz'); 
xlim([0 1e-3]); ylim([-2e4 6e4]);
klegend({'GRAPE (half time) $^{1}$H',...
         'GRAPE (half time) $^{15}$N'},...
         'Location','South');
subplot(2,2,4); plot(time_axis,real(fid_c)); drawnow();
klegend({'Adiabatic, tan ramp','GRAPE, same time',...
         'GRAPE, half time'},'Location','NorthEast');

end

