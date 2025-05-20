% Simulation of XiX/TPPM DNP field profile in the steady state 
% with electron Rabi frequency ensemble.
% 
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% i.kuprov@soton.ac.uk
% guinevere.mathies@uni-konstanz.de

function tppm_field_profile_ensemble_b1()

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

% B1 ensemble
[b1,wb1]=gaussleg(10e6,20e6,5); % Hz

% Microwave resonance offsets, Hz
offsets=linspace(-100e6,100e6,201);

% Experiment parameters
parameters.spins={'E','1H'};
parameters.grid='rep_2ang_800pts_sph';
parameters.pulse_dur=48e-9;              % Pulse duration, seconds
parameters.nloops=32;                    % Number of XiX/TPPM DNP blocks (power of 2)
parameters.phase=pi;                     % Second pulse inverted phase
parameters.shot_spacing=204e-6;
parameters.addshift=-13e6;
parameters.el_offs=offsets;

% Preallocate equilibrium DNP value array
dnp=zeros([numel(offsets) numel(b1)],'like',1i);

% Over B1 fields
for k=1:numel(b1)

    % Set electron nutation frequency
    parameters.irr_powers=b1(k);

    % Run the steady state simulation
    dnp(:,k)=powder(spin_system,@xixdnp_steady,parameters,'esr');

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);

% Plotting 
figure(); plot(parameters.el_offs/1e6,real(dnp)); 
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Microwave resonance offset, MHz');
kgrid; xlim tight;

end

