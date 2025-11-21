% Simulation of TOP DNP contact time dependence in the 
% steady state with electron Rabi frequency ensemble.
%
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function top_q_con_time_ensemble_b1()

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

% TOP loop count
loop_counts=1:256;

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
parameters.pulse_dur=10e-9;              % Pulse duration, seconds
parameters.delay_dur=14e-9;
parameters.addshift=-13e6;

% B1 ensemble A
[b1,wb1]=gaussleg(10e6,20e6,5); % Hz

% Preallocate equilibrium DNP value array
dnp_a=zeros([numel(loop_counts) numel(b1)],'like',1i);

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

        % Parameter set A
        localpar.el_offs=95e6;
        pulses_dur=localpar.nloops*(localpar.pulse_dur+...
                                    localpar.delay_dur);
        localpar.shot_spacing=102e-6 - pulses_dur;

        % Run the steady state simulation A
        dnp_a(m,k)=powder(spin_system,@topdnp_steady,localpar,'esr');

    end

end

% Integrate over the B1 field distribution
dnp_a=sum(dnp_a.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);

% B1 ensemble B
[b1,wb1]=gaussleg(25e6,35e6,5); % Hz

% Preallocate equilibrium DNP value array
dnp_b=zeros([numel(loop_counts) numel(b1)],'like',1i);

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

        % Parameter set B
        localpar.el_offs=92e6;
        pulses_dur=localpar.nloops*(localpar.pulse_dur+...
                                    localpar.delay_dur);
        localpar.shot_spacing=153e-6 - pulses_dur;

        % Run the steady state simulation B
        dnp_b(m,k)=powder(spin_system,@topdnp_steady,localpar,'esr');

    end

end

% Integrate over the B1 field distribution
dnp_b=sum(dnp_b.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);

% Plotting 
contact_times=(parameters.pulse_dur+parameters.delay_dur)*loop_counts;
kfigure(); plot(contact_times*1e6,real(dnp_a));
hold on; plot(contact_times*1e6,real(dnp_b));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');
klegend({'TOP, 18 MHz','TOP, 33 MHz'});
kxlabel('Total contact time, $\mu$s'); 
kgrid; xlim tight; ylim padded;

% Save the figure
savefig(gcf,'top_q_con_time_ensemble_b1.fig');

end

