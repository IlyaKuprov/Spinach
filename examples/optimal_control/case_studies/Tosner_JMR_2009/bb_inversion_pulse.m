% Broadband inversion pulse design for liquid-state NMR. Reprodu-
% ces, using Spinach, the second example from:
%
%             http://dx.doi.org/10.1016/j.jmr.2008.11.020
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
Sz=state(spin_system,'Lz',1);
Sz=Sz/norm(full(Sz),2);

% Control and offset operators
LxH=operator(spin_system,'Lx',1);
LyH=operator(spin_system,'Ly',1);
LzH=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Control data structure
control.drifts={{H}};                            % Drift Hamiltonian
control.operators={LxH,LyH};                     % Control operators
control.off_ops={LzH};                           % Offset operator
control.offsets={linspace(-50e3,50e3,101)};      % As per the paper
control.rho_init={+Sz};                          % Initial state
control.rho_targ={-Sz};                          % Target state
control.pulse_dt=1e-6*ones(1,600);               % As per the paper
control.pwr_levels=2*pi*10e3;                    % As per the paper
control.penalties={'NS','SNSA'};                 % Penalties
control.p_weights=[0.01 10];                     % Penalty weights
control.method='lbfgs';                          % Optimiser
control.max_iter=200;                            % Max iterations
control.parallel='ensemble';                     % Parallel mode
control.plotting={'phi_controls','amp_controls',...
                  'robustness','spectrogram'};

% Random guess
guess=randn(2,600)/10;

% Optimisation
spin_system=optimcon(spin_system,control);
xy_profile=fmaxnewton(spin_system,@grape_xy,guess);

% Return to physical units
rf_scale=mean(control.pwr_levels);
CLx=rf_scale*xy_profile(1,:);
CLy=rf_scale*xy_profile(2,:);

% Offset grid for verification
offs_hz=linspace(-100e3,100e3,201);
inv_eff=zeros(size(offs_hz));

% Test simulation
parfor k=1:numel(offs_hz)
    
    % Add the offset term
    Hd=H+2*pi*offs_hz(k)*LzH; 

    % Run the pulse
    rho_f=shaped_pulse_xy(spin_system,Hd,{LxH, LyH},{CLx,CLy}, ...
                          control.pulse_dt,Sz,'expv-pwc');     %#ok<PFBNS>

    % Compute inversion efficiency
    inv_eff(k)=-real(Sz'*rho_f);

end

% Plot the efficiency profile
kfigure(); plot(offs_hz/1e3,inv_eff); kgrid;
kxlabel('offset, kHz'); kylabel('fidelity');
xlim([-100 100]); ylim([-1.1 1.1]);

end

