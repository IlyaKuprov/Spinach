% Simulation of XiX DNP field profile in the steady state,
% a single spin system without ensemble averaging.
% 
% Calculation time: seconds
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function xix_w_field_profile_single()

% W-band magnet
sys.magnet=3.4;

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
xyz=cell2mat(inter.coordinates); r_en=xyz(2,3);

% Relaxation rates, distance and ori. dep. R1n
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1e-3,52,r_en,bet); 
inter.r1_rates={1e3 r1n_rate};
inter.r2_rates={200e3 50e3};
inter.rlx_keep='diagonal';
inter.equilibrium='dibari';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'hygiene'}';
sys.enable={'op_cache','ham_cache'};

% Propagator accuracy
sys.tols.prop_chop=1e-12;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Detect the proton
parameters.coil=state(spin_system,'Lz',2);

% Experiment parameters
parameters.spins={'E','1H'};
parameters.irr_powers=20e6;              % Electron nutation frequency [Hz]
parameters.pulse_dur=18e-9;              % Pulse duration, seconds
parameters.grid='rep_2ang_800pts_sph';
parameters.nloops=10;                    % Number of XiX DNP blocks
parameters.phase=pi;                     % Second pulse inverted phase
parameters.shot_spacing=167e-6;
parameters.addshift=-33e6;
parameters.el_offs=linspace(-300e6,300e6,201);

% Run the steady state simulation
dnp=powder(spin_system,@xixdnp_steady,parameters,'esr');

% Plotting 
figure(); plot(parameters.el_offs/1e6,real(dnp)); 
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Microwave resonance offset, MHz'); 
kgrid; xlim tight; ylim padded;

end

