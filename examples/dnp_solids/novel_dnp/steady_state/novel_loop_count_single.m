% Simulation of NOVEL DNP loop count dependence in the 
% steady state.
%
% Calculation time: seconds.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function novel_loop_count_single()

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
inter.r1_rates={1000 r1n_rate};
inter.r2_rates={200000 50e3};
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

% Number of NOVEL loops
loop_counts=[1 2 4 8 16 32 64 128 256];

% Experiment parameters
parameters.spins={'E','1H'};
parameters.irr_powers=17.8e6;            % Electron nutation frequency [Hz]
parameters.grid='rep_2ang_800pts_sph';
parameters.contact_dur=500e-9;           % Pulse duration, seconds
parameters.flippulse=1;                  % 1 for NOVEL, 0 for SE
parameters.flipback=1;                   % 1 for flipback, 0 for no flipback
parameters.addshift=-3.3e6;
parameters.el_offs=0e6;
       
% Over loop counts
dnp=zeros(size(loop_counts),'like',1i);
for m=1:numel(loop_counts)

    % Set the number of loops
    parameters.nloops=loop_counts(m);

    % Run the steady state simulation
    dnp(m)=powder(spin_system,@noveldnp_steady,parameters,'esr');

end

% Plotting 
contact_times=parameters.pulse_dur*2*loop_counts;
figure(); plot(contact_times*1e6,real(dnp),'-o');
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Contact time, $\mu$s'); kgrid; xlim tight;

end

