% Simulation of XiX DNP field profile in the steady state 
% with electron-proton distance ensemble averaging.
% 
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function xix_q_field_profile_ensemble_r()

% Q-band magnet
sys.magnet=1.2142;

% Electron and proton
sys.isotopes={'E','1H'};

% Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]};
inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]};

% Spin temperature
inter.temperature=80;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Propagator accuracy
sys.tols.prop_chop=1e-12;

% Algorithmic options
sys.disable={'hygiene'}';
sys.enable={'op_cache','ham_cache'};

% Distance ensemble, Gauss-Legendre points
[r,w]=gaussleg(3.5,20,3);

% Microwave resonance offsets, Hz
offsets=linspace(-100e6,100e6,201);

% Compute DNP at each distance
dnp=zeros([numel(offsets) numel(r)]);
for n=1:numel(r)

    % Cartesian coordinates
    inter.coordinates={[0.000 0.000 0.000];
                       [0.000 0.000 r(n) ]};
    
    % Relaxation rates, distance and ori. dep. R1n
    inter.relaxation={'t1_t2'};
    r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                                   2.00230,1e-3,52,r(n),bet); 
    inter.r1_rates={1e3 r1n_rate};
    inter.r2_rates={200e3 50e3};
    inter.rlx_keep='diagonal';
    inter.equilibrium='dibari';

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Detect the proton
    parameters.coil=state(spin_system,'Lz',2);

    % Experiment parameters
    parameters.spins={'E','1H'};
    parameters.irr_powers=18e6;              % Electron nutation frequency [Hz]
    parameters.grid='rep_2ang_800pts_sph';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=36;                    % Number of XiX DNP blocks (power of 2)
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.shot_spacing=204e-6;
    parameters.addshift=-13e6;
    parameters.el_offs=offsets;

    % Run the steady state simulation
    dnp(:,n)=powder(spin_system,@xixdnp_steady,parameters,'esr');

end

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(w,[1 numel(w)]),2)/sum((r.^2).*w);
        
% Plotting 
figure(); plot(parameters.el_offs/1e6,real(dnp)); 
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Microwave resonance offset, MHz'); 
kgrid; xlim tight; ylim padded;

end

