% Spinach variant of the broadband refocusing example from
%
%          http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Unlike the paper, this script uses a 30 kHz nominal RF scale
% with a soft Cartesian penalty, and scores three state transfers
% rather than the published target propagator.
%
% GRAPE is used to design a 200 µs broadband x-phase π pulse:
%
%              {Sx ->  Sx,  Sy -> -Sy,  Sz -> -Sz}
%
% over an offset range of ±12.5 kHz.
%
% aditya.dev@weizmann.ac.il

function bb_refocusing_pulse()

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

% Normalised Cartesian basis states
rho_x=state(spin_system,'Lx',1); rho_x=rho_x/norm(full(rho_x),2);
rho_y=state(spin_system,'Ly',1); rho_y=rho_y/norm(full(rho_y),2);
rho_z=state(spin_system,'Lz',1); rho_z=rho_z/norm(full(rho_z),2);

% RF controls and offset operator
lx=operator(spin_system,'Lx',1);
ly=operator(spin_system,'Ly',1);
lz=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Control parameters
control.isotopes={'1H'};                         % Isotopes
control.channels=[1; 1];                         % Channel map
control.drifts={{H}};                            % Drift Hamiltonian
control.operators={lx,ly};                       % Control operators
control.off_ops={lz};                            % Offset operator
control.offsets={linspace(-12.5e3,12.5e3,101)};  % As per the paper, Hz
control.rho_init={rho_x,rho_y,rho_z};            % Initial states
control.rho_targ={rho_x,-rho_y,-rho_z};          % Target states
control.pulse_dt=(200e-6/600)*ones(1,600);       % 200 us, 600 Spinach slices
control.pwr_levels=2*pi*30e3;                    % Spinach scale; paper says 15 kHz
control.penalties={'NS','SNS'};                  % Penalties
control.p_weights=[0.01 100];                    % Penalty weights
control.method='lbfgs';                          % Optimiser
control.max_iter=200;                            % Max iterations

% Visual diagnostics
control.plotting={'phi_controls','xy_controls',...
                  'spectrogram','robustness'};

% Random initial guess
guess=randn(2,600)/10;

% Spinach optimal-control housekeeping
spin_system=optimcon(spin_system,control);

% Run the optimisation
xy_profile=fmaxnewton(spin_system,@grape_xy,guess);

% Convert normalised waveform to physical rad/s controls
rf_x=control.pwr_levels*xy_profile(1,:);
rf_y=control.pwr_levels*xy_profile(2,:);

% Offset grid for verification
offs_hz=linspace(-25e3,25e3,201);
fidelities=zeros(size(offs_hz));

% Test simulation
parfor k=1:numel(offs_hz)

    % Get drift Hamiltonian
    drift=H+2*pi*offs_hz(k)*lz;

    % Apply the pulse
    rho_final=shaped_pulse_xy(spin_system,drift,control.operators,{rf_x,rf_y},...
                              control.pulse_dt,[rho_x,rho_y,rho_z],'expv-pwc'); %#ok<PFBNS>

    % Calculate the fidelity
    fidelities(k)=real(trace([rho_x,-rho_y,-rho_z]'*rho_final))/3;

end

% Plot the fidelity profile
kfigure(); plot(offs_hz/1e3,fidelities); kgrid;
kxlabel('offset, kHz'); kylabel('fidelity');
xlim([-25 25]); ylim([-1.1 1.1]);

end


