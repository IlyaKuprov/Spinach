% Simulation of TOP DNP contact time dependence in the 
% steady state.
%
% Calculation time: seconds.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function top_q_con_time_single()

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

% Propagator accuracy
sys.tols.prop_chop=1e-12;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Detect the proton
parameters.coil=state(spin_system,'Lz',2);

% TOP loop count
loop_counts=1:256;

% Experiment parameters
parameters.spins={'E','1H'};
parameters.grid='rep_2ang_800pts_sph';
parameters.pulse_dur=10e-9;              % Pulse duration, seconds
parameters.delay_dur=14e-9;
parameters.addshift=-13e6;
       
% Over loop counts
dnp_a=zeros(size(loop_counts),'like',1i);
dnp_b=zeros(size(loop_counts),'like',1i);
parfor m=1:numel(loop_counts)

    % Localise parameters
    localpar=parameters;

    % Set the number of loops
    localpar.nloops=loop_counts(m);

    % Parameter set A
    localpar.irr_powers=18e6;             
    localpar.el_offs=95e6;
    pulses_dur=localpar.nloops*(localpar.pulse_dur+...
                                localpar.delay_dur);
    localpar.shot_spacing=102e-6 - pulses_dur;
    
    % Run the steady state simulation A
    dnp_a(m)=powder(spin_system,@topdnp_steady,localpar,'esr');

    % Parameter set B
    localpar.irr_powers=33e6;             
    localpar.el_offs=92e6;
    pulses_dur=localpar.nloops*(localpar.pulse_dur+...
                                localpar.delay_dur);
    localpar.shot_spacing=153e-6 - pulses_dur;
    
    % Run the steady state simulation B
    dnp_b(m)=powder(spin_system,@topdnp_steady,localpar,'esr');

end

% Plotting 
contact_times=(parameters.pulse_dur+parameters.delay_dur)*loop_counts;
figure(); plot(contact_times*1e6,real(dnp_a));
hold on; plot(contact_times*1e6,real(dnp_b));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');
klegend({'TOP, 18 MHz','TOP, 33 MHz'});
kxlabel('Total contact time, $\mu$s'); 
kgrid; xlim tight; ylim padded;

% Save for later
savefig(gcf,'top_q_con_time_single.fig');
save('top_q_con_time_single.mat');

end

