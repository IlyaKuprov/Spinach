% Simulation of XiX DNP repetition time scan in the 
% steady state.
% 
% Calculation time: seconds.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function xix_rep_time_single()

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
sys.enable={'op_cache','ham_cache'};

% Shot spacings, s
srt=logspace(-5,-3,30);

% Preallocate equilibrium DNP value array
dnp=zeros(size(srt),'like',1i);

% Relaxation rates, distance and orientation
% dependence provided using a function handle
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1e-3,52,r_en,bet);
inter.r1_rates={1000 r1n_rate};
inter.r2_rates={200000 50e3};
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
parameters.irr_powers=17.8e6;            % Electron nutation frequency [Hz]
parameters.pulse_dur=48e-9;              % Pulse duration, seconds
parameters.nloops=32;                    % Number of XiX DNP blocks (power of 2)
parameters.phase=pi;                     % Second pulse inverted phase
parameters.addshift=-13e6;
parameters.el_offs=-39e6;

% Over shot spacing
for m=1:numel(srt)

    % Set the shot spacing
    parameters.shot_spacing=srt(m);

    % Run the steady state simulation
    dnp(m)=powder(spin_system,@xixdnp_steady,parameters,'esr');

end

% Plotting 
figure(); plot(srt*1e3,real(dnp));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Repetition time, ms');
kgrid; xlim tight;

end

