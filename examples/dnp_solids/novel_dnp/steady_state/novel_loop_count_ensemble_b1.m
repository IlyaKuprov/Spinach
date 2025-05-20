% Simulation of NOVEL DNP loop count dependence in the 
% steady state with electron Rabi frequency ensemble.
%
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function novel_loop_count_ensemble_b1()

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

% Distance and B1 ensemble
[b1,wb1]=gaussleg(14e6,16e6,5); % Hz

% Number of NOVEL loops
loop_counts=[1 2 4 8 16 32 64 128 256];

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
parameters.contact_dur=500e-9;           % Pulse duration, seconds
parameters.flippulse=1;                  % 1 for NOVEL, 0 for SE
parameters.flipback=1;                   % 1 for flipback, 0 for no flipback
parameters.addshift=-3.3e6;
parameters.el_offs=0e6;

% Preallocate equilibrium DNP value array
dnp=zeros([numel(loop_counts) numel(b1)],'like',1i);

% Over B1 fields
for k=1:numel(b1)

    % Set electron nutation frequency
    parameters.irr_powers=b1(k);

    % Over loop counts
    for m=1:numel(loop_counts)

        % Set the number of loops
        parameters.nloops=loop_counts(m);

        % Run the steady state simulation
        dnp(m,k)=powder(spin_system,@noveldnp_steady,parameters,'esr');

    end

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);

% Plotting 
contact_times=parameters.pulse_dur*2*loop_counts;
figure(); plot(contact_times*1e6,real(dnp),'-o');
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Contact time, $\mu$s'); kgrid; xlim tight;

end

