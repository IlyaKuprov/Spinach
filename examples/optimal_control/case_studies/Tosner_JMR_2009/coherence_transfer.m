% The first optimal control example from 
%
%       http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% A heteronuclear two-spin system (1H–13C) with an scalar
% and both nuclei set on resonance; the goal is to trans-
% fer transverse magnetisation from proton to carbon:
%
%                         Hx → Cx
%
% over a fixed evolution period T = 1/J.
%
% aditya.dev@weizmann.ac.il

function coherence_transfer()

% Magnetic field, Tesla
sys.magnet=14.1;
sys.isotopes={'1H','13C'};

% Chemical shifts, ppm
inter.zeeman.scalar={0.0,0.0};

% Scalar coupling, Hz
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=140;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system, bas);

% Initial state: Lx on proton (spin 1)
rho_init=state(spin_system,{'Lx'},{1});
rho_init=rho_init/norm(rho_init,2);

% Target state: Lx on carbon (spin 2)
rho_targ=state(spin_system,{'Lx'},{2});
rho_targ=rho_targ/norm(full(rho_targ), 2);

% Control operators
LxH=operator(spin_system,'Lx',1);
LyH=operator(spin_system,'Ly',1);
LxC=operator(spin_system,'Lx',2);
LyC=operator(spin_system,'Ly',2);

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Control data structure
control.drifts={{H}};                          % Drift Hamiltonian
control.operators={LxH,LyH,LxC,LyC};           % Control operators
control.rho_init={rho_init};                   % Initial state
control.rho_targ={rho_targ};                   % Target state
control.pulse_dt=(1/(140*150))*ones(1,150);    % Duration 1/J, 150 points
control.pwr_levels=2*pi*linspace(10,1000,10);  % As per the paper
control.penalties={'NS'};                      % Penalties
control.p_weights=0.01;                        % Penalty weight
control.method='lbfgs';                        % Optimiser
control.max_iter=200;                          % Max iterations
control.parallel='ensemble';                   % Parallel mode

% Visual diagnostics
control.plotting={'xy_controls','spectrogram','robustness'};

% Random initial guess
guess=rand(4,150)/10; 

% Optimisation
spin_system=optimcon(spin_system,control);
xy_profile=fmaxnewton(spin_system,@grape_xy,guess);

% Better plotting needed here

end
