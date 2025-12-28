% Broadband inversion pulse design for liquid-state NMR. This script
% reproduces,using Spinach,the second SIMPSON optimal control exam-
% ple from http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% A single-spin (1H) system is considered in the rotating frame,on
% resonance (no isotropic chemical shift). The goal is to design a
% broadband inversion pulse that performs:
%
%                              I_z  →  -I_z
%
% uniformly over a frequency offset range of ±50 kHz,under a fixed 
% pulse duration T=600 µs,discretised into 600 time steps of 1 µs.
%
% The control fields are Cartesian RF components (Lx,Ly) on 1H,and
% robustness is enforced by an offset ensemble and RF-power penalties.
%
% aditya.dev@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function broadband_inversion_pulse_for_nmr()

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
control.offsets={2*pi*linspace(-50e3,50e3,101)}; % As per the paper
control.rho_init={Sz};                           % Initial state
control.rho_targ={-Sz};                          % Target state
control.pulse_dt=1e-6*ones(1,600);               % As per the paper
control.pwr_levels=2*pi*13e3;                    % As per the paper
control.penalties={'NS','SNS'};                  % Penalties
control.p_weights=[0.1 100];                     % Penalty weights
control.method='lbfgs';                          % Optimiser
control.max_iter=200;                            % Max iterations
control.parallel='ensemble';                     % Parallel mode
control.plotting={'xy_controls','amp_controls','robustness'};

% Random guess
guess=randn(2,600)/3;

% Optimisation
spin_system=optimcon(spin_system,control);
xy_profile=fminnewton(spin_system,@grape_xy,guess);

% Demonstration simulation goes here

end

