% Simulation of TOP DNP repetition time scan in the 
% steady state.
% 
% Calculation time: seconds.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function top_q_rep_time_single()

% Q-band magnet
sys.magnet=1.2142;

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

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Propagator accuracy
sys.tols.prop_chop=1e-12;

% Algorithmic options
sys.disable={'hygiene'}';

% Relaxation rates, distance and orientation
% dependence provided using a function handle
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1e-3,52,r_en,bet);
inter.r1_rates={1e3 r1n_rate};
inter.r2_rates={200e3 50e3};
inter.rlx_keep='diagonal';
inter.equilibrium='dibari';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Detect the proton
parameters.coil=state(spin_system,'Lz','1H');

% Experiment parameters
parameters.spins={'E','1H'};
parameters.grid='rep_2ang_800pts_sph';
parameters.irr_powers=18e6;              % Electron nutation frequency [Hz]
parameters.pulse_dur=10e-9;              % Pulse duration, seconds
parameters.delay_dur=14e-9;              % Delay duration, seconds
parameters.nloops=300;                   % Number of TOP DNP blocks
parameters.addshift=-13e6;
parameters.el_offs=95e6;

% Log spacing for rep. time
rep_time=logspace(-5,-3,30);

% Preallocate equilibrium DNP value array
dnp=zeros(size(rep_time),'like',1i);

% Over repetition times
parfor m=1:numel(rep_time)

    % Localise parameters
    localpar=parameters;

    % Set the shot spacing
    pulses_time=localpar.nloops*(localpar.pulse_dur+...
                                 localpar.delay_dur);
    localpar.shot_spacing=rep_time(m)-pulses_time;

    % Run the steady state simulation
    dnp(m)=powder(spin_system,@topdnp_steady,localpar,'esr');

end

% Plotting 
figure(); plot(rep_time*1e3,real(dnp));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Repetition time, ms'); 
kgrid; xlim([0 2]); ylim padded;

% Save for later
savefig(gcf,'top_q_rep_time_single.fig');

end

