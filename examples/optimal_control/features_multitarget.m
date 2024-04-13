% An example of multi-target optimal control pulse design in the context
% of singlet state NMR spectroscopy. A pulse is designed that moves TT 
% (carbon-triplet, proton-triplet) into SS (carbon-singlet, proton-singlet)
% and TS (carbon-triplet, proton-singlet) into ST (carbon-singlet, proton-
% triplet) simultaneously. The system is assumed to have a distribution in
% one of the J-couplings.
%
% Calculation time: hours.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function features_multitarget()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','13C','13C','1H'};

% Interactions
inter.zeeman.scalar={0.0, 0.0, 0.0, 0.0};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=15.0;
inter.coupling.scalar{3,4}=15.0;
inter.coupling.scalar{1,3}=3.0;
inter.coupling.scalar{2,4}=3.0;
inter.coupling.scalar{2,3}=150;
inter.coupling.scalar{1,4}=8.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none'; 

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up relevant operators
unit_left=unit_oper(spin_system);
HxHx_left=operator(spin_system,{'Lx','Lx'},{1,4},'left');
HyHy_left=operator(spin_system,{'Ly','Ly'},{1,4},'left');
HzHz_left=operator(spin_system,{'Lz','Lz'},{1,4},'left');
CxCx_left=operator(spin_system,{'Lx','Lx'},{2,3},'left');
CyCy_left=operator(spin_system,{'Ly','Ly'},{2,3},'left');
CzCz_left=operator(spin_system,{'Lz','Lz'},{2,3},'left');
HS_left=unit_left/4-(HzHz_left+HyHy_left+HxHx_left);
CS_left=unit_left/4-(CzCz_left+CyCy_left+CxCx_left);
HT_left=unit_left/4-(HzHz_left-HyHy_left-HxHx_left);
CT_left=unit_left/4-(CzCz_left-CyCy_left-CxCx_left);

% Set up source and target states
TT=CT_left*HT_left*unit_state(spin_system); TT=TT/norm(TT,2);
TS=CT_left*HS_left*unit_state(spin_system); TS=TS/norm(TS,2);
ST=CS_left*HT_left*unit_state(spin_system); ST=ST/norm(ST,2);
SS=CS_left*HS_left*unit_state(spin_system); SS=SS/norm(SS,2);
sources={TT,TS}; targets={SS,ST};

% Get the control operators
LxH=operator(spin_system,'Lx','1H');
LyH=operator(spin_system,'Ly','1H');
LxC=operator(spin_system,'Lx','13C'); 
LyC=operator(spin_system,'Ly','13C');

% Generate a distribution in J-couplings
J=linspace(13,16,11); H=cell(numel(J),1);
for m=1:numel(H)
    spin_system.inter.coupling.matrix{1,2}=2*pi*J(m)*eye(3);
    H{m}{1}=hamiltonian(assume(spin_system,'nmr'));
end

% Define control parameters
control.drifts=H;                       % Drift array
control.operators={LxH,LyH,LxC,LyC};	% Control operators
control.rho_init=sources;               % Initial states
control.rho_targ=targets;               % Destination states
control.pwr_levels=2*pi*500;            % Power levals
control.pulse_dt=5e-4*ones(1,275);      % Interval grid
control.penalties={'SNS'};              % Penalty
control.p_weights=100;                  % Penalty weight
control.method='lbfgs';                 % Optimisation method
control.max_iter=150;                   % Termination condition
control.parallel='ensemble';            % Parallelisation

% Control trajectory analysis plots
control.plotting={'correlation_order','local_each_spin',...
                  'xy_controls','robustness','spectrogram'};

% Make control system structure
spin_system=optimcon(spin_system,control);

% Run the optimisation from a random guess
pulse=fminnewton(spin_system,@grape_xy,randn(4,275)/3);

% Denormalise and format the waveform
pulse=pulse*control.pwr_levels;
pulse=mat2cell(pulse,[1 1 1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,H{6}{1},control.operators,pulse,...
                    control.pulse_dt,cell2mat(control.rho_init),'expv-pwc');
fidelity_a=real(control.rho_targ{1}'*rho(:,1));
fidelity_b=real(control.rho_targ{2}'*rho(:,2));
report(spin_system,['TT->SS: Re[<target|rho(T)>] = ' num2str(fidelity_a)]);
report(spin_system,['TS->ST: Re[<target|rho(T)>] = ' num2str(fidelity_b)]);

end

