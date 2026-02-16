% Spinach implementation of the broadband refocusing example from
%
%          http://dx.doi.org/10.1016/j.jmr.2008.11.020
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
Sx=state(spin_system,'Lx',1); Sx=Sx/norm(full(Sx),2);
Sy=state(spin_system,'Ly',1); Sy=Sy/norm(full(Sy),2);
Sz=state(spin_system,'Lz',1); Sz=Sz/norm(full(Sz),2);

% RF controls and offset operator
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);
Lz=operator(spin_system,'Lz',1);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Control data structure
control.drifts={{H}};                           % Drift Hamiltonian
control.operators={Lx,Ly};                      % Control operators
control.off_ops={Lz};                           % Offset operator
control.offsets={linspace(-12.5e3,12.5e3,101)}; % As per the paper, Hz
control.rho_init={Sx, Sy, Sz};                  % Initial states
control.rho_targ={Sx,-Sy,-Sz};                  % Target states
control.pulse_dt=(200e-6/600)*ones(1,600);      % As per the paper
control.pwr_levels=2*pi*30e3;                   % As per the paper
control.penalties={'NS','SNS'};                 % Penalties
control.p_weights=[0.01 100];                   % Penalty weights
control.method='lbfgs';                         % Optimiser
control.max_iter=200;                           % Max iterations
control.parallel='ensemble';                    % Parallel mode

% Visual diagnostics
control.plotting={'phi_controls','xy_controls',...
                  'spectrogram','robustness'};

% Random initial guess
guess=randn(2,600)/10;

% Optimisation
spin_system=optimcon(spin_system,control);
xy_profile=fmaxnewton(spin_system,@grape_xy,guess);

% Convert normalised waveform to physical rad/s controls
CLx=control.pwr_levels*xy_profile(1,:);
CLy=control.pwr_levels*xy_profile(2,:);

% Test the pulse 
offs_hz=linspace(-25e3,25e3,201);
fidelities=zeros(size(offs_hz));
parfor k=1:numel(offs_hz)

    % Get drift Hamiltonian
    Hk=H+2*pi*offs_hz(k)*Lz;

    % Apply the pulse
    rho_k=shaped_pulse_xy(spin_system,Hk,control.operators,{CLx,CLy}, ...
                          control.pulse_dt,[Sx,Sy,Sz],'expv-pwc'); %#ok<PFBNS>

    % Calculate the fidelity
    fidelities(k)=real(trace([Sx,-Sy,-Sz]'*rho_k))/3;

end

% Plot the fidelity profile
kfigure(); plot(offs_hz/1e3,fidelities); kgrid;
kxlabel('offset, kHz'); kylabel('fidelity');
xlim([-25 25]); ylim([-1.1 1.1]);

end

