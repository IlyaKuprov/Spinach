% Simulation of TPPM DNP contact time dependence in the 
% steady state with electron Rabi frequency ensemble.
%
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function tppm_q_con_time_ensemble_b1()

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

% Distance and B1 ensemble
[b1,wb1]=gaussleg(25e6,35e6,5); % Hz

% TPPM loop count
loop_counts=1:256;

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
parameters.irr_powers=33e6;              % Electron nutation frequency [Hz]
parameters.pulse_dur=16e-9;              % Pulse duration, seconds
parameters.grid='rep_2ang_800pts_sph';
parameters.phase=120*pi/180;             % Second pulse phase
parameters.addshift=-13e6;
parameters.el_offs=2e6;

% Preallocate equilibrium DNP value array
dnp=zeros([numel(loop_counts) numel(b1)],'like',1i);

% Over B1 fields
for k=1:numel(b1)

    % Set electron nutation frequency
    parameters.irr_powers=b1(k);

    % Over loop counts
    parfor m=1:numel(loop_counts)

        % Localise parameters
        localpar=parameters;

        % Set the number of loops
        localpar.nloops=loop_counts(m);

        % Update the shot spacing
        pulses_dur=2*localpar.nloops*localpar.pulse_dur;
        localpar.shot_spacing=816e-6 - pulses_dur;

        % Run the steady state simulation
        dnp(m,k)=powder(spin_system,@xixdnp_steady,localpar,'esr');

    end

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);

% Plotting 
contact_times=parameters.pulse_dur*2*loop_counts;
figure(); plot(contact_times*1e6,real(dnp));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Total contact time, $\mu$s'); 
kgrid; xlim tight; ylim padded;

% Save for later
savefig(gcf,'tppm_q_con_time_ensemble_b1.fig');
save('tppm_q_con_time_ensemble_b1.mat');

end

