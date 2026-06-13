% Panoramic optimisation for stroboscopic steady state DNP
% with the timing and power settings matching the XiX ex-
% eriment, but complete liberty is the choice of phase.
%
% guinevere.mathies@uni-konstanz.de
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function solid_effect_xix()

% W-band magnet
sys.magnet=3.35316; % HIPER at St Andrews

% Electron and proton
sys.isotopes={'E','1H'};

% Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]};
inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]};

% Spin temperature
inter.temperature=80;

% Cartesian coordinates
inter.coordinates={[0.000 0.000 0.000];
                   [0.000 0.000 3.500]};

% Get electron-nuclear distance
r_en=norm(inter.coordinates{2}-...
          inter.coordinates{1},2);

% Relaxation rates, distance and ori. dep. R1n
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1.0e-3,52.0,r_en,bet); 
inter.r1_rates={1e3 r1n_rate};
inter.r2_rates={200e3 50e3};
inter.rlx_keep='diagonal';
inter.equilibrium='dibari';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Parallelisation settings
sys.parallel={'processes',240};

% Very tight numerics
sys.tols.prop_chop=1e-14;
sys.tols.stst_tol=1e-10;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get electron control operators
Ex=operator(spin_system,'Lx','E');
Ey=operator(spin_system,'Ly','E');

% Get electron offset operator
Ez=operator(spin_system,'Lz','E');

% Initial state will be ignored by StStSt, but set to
% thermodynamic equilibrium for when the option is off
H=hamiltonian(assume(spin_system,'labframe'),'left');
rho_init=equilibrium(spin_system,H);

% Target state is nuclear magnetisation,
% relative to the thermal equilibrium
rho_targ=state(spin_system,'Lz','1H');
thermal=real(rho_targ'*rho_init);
rho_targ=rho_targ/thermal;

% Powder context parameters
parameters.spins={'E','1H'};
parameters.grid='rep_2ang_800pts_sph';

% Transmitter is at 94.0 GHz precisely
parameters.offset=[-spin('E')*sys.magnet/(2*pi)-94.0e9, 0];

% Get drift Liouvillians at all orientations
control.drifts=drifts(spin_system,@powder,parameters,'esr');

% Define control parameters
control.operators={Ex,Ey};                       % Controls
control.rho_init={rho_init};                     % Starting state
control.rho_targ={rho_targ};                     % Destination state
control.pwr_levels=2*pi*linspace(5,25,20)*1e6;   % Microwave power, rad/s
control.off_ops={Ez};                            % Offset operator
control.offsets={1e6*[-2 -1  0 +1 +2]};          % Offsets, Hz
control.pulse_dt=[0.5e-9*ones(1,720) ...         % Pulse itself
                  0.5e-9*ones(1,20)  ...         % Ringdown delay
                  167e-6];                       % Sequence delay
control.freeze=[false(1,720) ...                 % Pulse itself
                true(1,20)   ...                 % Ringdown delay
                true(1,1)];                      % Sequence delay
control.amplitudes=[ones(1,720) ...              % Pulse itself
                    zeros(1,20) ...              % Ringdown delay
                    zeros(1,1)];                 % Sequence delay
control.method='rbfgs';                          % Optimisation method
control.max_iter=10000;                          % Maximum iterations
control.parallel='ensemble';                     % Parallelisation
control.steady=true();                           % Steady state
control.budget=500;                              % Increase this

% Plotting options
control.plotting={'robustness','spectrogram'};

% Load HiPER filter function
load('hiper_kernel_trans.mat','h'); 
h=h(1:16); h=h/abs(sum(h));                      % Unit abs DC gain
control.distortion={@(w)firf(w,h)};              % For optimisation 
control.distplot={@(w)firf(w,h)};                % For plotting

% Optimal control housekeeping
spin_system=optimcon(spin_system,control);

% XiX waveform, shifted up by 160 MHz
altern=72*(1:1:10);
for n=1:numel(altern)
    j=altern(n); guess(j-37:j)=pi;
end
guess=guess-wrapTo2Pi(2*pi*140e6*linspace(0,360e-9,720));
guess=[guess zeros(1,20) 0];

% Run the optimisation and save the entire workspace
pulse_profile=fmaxnewton(spin_system,@grape_phase,guess); %#ok<NASGU>

end

