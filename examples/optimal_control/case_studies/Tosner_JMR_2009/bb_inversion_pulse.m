% Spinach variant of the broadband inversion example from:
%
%             http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% SNSA penalises amplitudes above 10 kHz rather than imposing the
% hard RF cap used in the article; this is not its published pulse.
%
% A single proton is considered in the rotating frame with multiple
% transmitter offsets (or chemical shifts). The goal is to design a
% 600 µs broadband inversion pulse (1 µs slices) that performs:
%
%                            I_z  →  -I_z
%
% uniformly over a frequency offset range of ±50 kHz; controls are
% Cartesian (Lx, Ly) operators in the rotating frame.
%
% aditya.dev@weizmann.ac.il

function bb_inversion_pulse()

% Magnetic field (Tesla)
sys.magnet=14.1;

% Isotope
sys.isotopes={'1H'};

% Chemical shift (ppm)
inter.zeeman.scalar={0.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial and target states
rho_z=state(spin_system,'Lz',1);
rho_z=rho_z/norm(full(rho_z),2);

% Control and offset operators
lx_h=operator(spin_system,'Lx',1);
ly_h=operator(spin_system,'Ly',1);
lz_h=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Control parameters
control.isotopes={'1H'};                     % Isotopes
control.channels=[1; 1];                     % Channel map
control.drifts={{H}};                        % Drift Hamiltonian
control.operators={lx_h,ly_h};               % Control operators
control.off_ops={lz_h};                      % Offset operator
control.offsets={linspace(-50e3,50e3,101)};  % As per the paper
control.rho_init={+rho_z};                   % Initial state
control.rho_targ={-rho_z};                   % Target state
control.pulse_dt=1e-6*ones(1,600);           % As per the paper
control.pwr_levels=2*pi*10e3;                % RF scale, not a hard cap
control.penalties={'NS','SNSA'};             % Penalties
control.p_weights=[0.01 10];                 % Penalty weights
control.method='lbfgs';                      % Optimiser
control.max_iter=200;                        % Max iterations
control.plotting={'phi_controls','amp_controls',...
                  'robustness','spectrogram'};

% Random guess
guess=randn(2,600)/10;

% Spinach optimal-control housekeeping
spin_system=optimcon(spin_system,control);

% Run the optimisation
xy_profile=fmaxnewton(spin_system,@grape_xy,guess);

% Return to physical units
rf_scale=mean(control.pwr_levels);
rf_x=rf_scale*xy_profile(1,:);
rf_y=rf_scale*xy_profile(2,:);

% Offset grid for verification
offs_hz=linspace(-100e3,100e3,201);
inv_eff=zeros(size(offs_hz));

% Test simulation
parfor k=1:numel(offs_hz)

    % Add the offset term
    drift=H+2*pi*offs_hz(k)*lz_h;

    % Run the pulse
    rho_final=shaped_pulse_xy(spin_system,drift,{lx_h,ly_h},{rf_x,rf_y},...
                              control.pulse_dt,rho_z,'expv-pwc');     %#ok<PFBNS>

    % Compute inversion efficiency
    inv_eff(k)=-real(rho_z'*rho_final);

end

% Plot the fidelity profile
kfigure(); plot(offs_hz/1e3,inv_eff); kgrid;
kxlabel('offset, kHz'); kylabel('fidelity');
xlim([-100 100]); ylim([-1.1 1.1]);

end


