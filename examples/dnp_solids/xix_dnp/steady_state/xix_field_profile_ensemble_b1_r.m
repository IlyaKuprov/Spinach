% Simulation of XiX DNP field profile in the steady state with
% electron-proton distance and electron Rabi frequency ensembles.
% 
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function xix_field_profile_ensemble_b1_r()

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

% Distance and B1 ensemble
[r,wr]=gaussleg(3.5,20,3);      % Angstrom
[b1,wb1]=gaussleg(10e6,20e6,5); % Hz

% Microwave resonance offsets, Hz
offsets=linspace(-100e6,100e6,201);

% Preallocate equilibrium DNP value array
dnp=zeros([numel(offsets) numel(r) numel(b1)],'like',1i);

% Over distances
for n=1:numel(r)

    % Cartesian coordinates
    inter.coordinates={[0.000 0.000 0.000];
                       [0.000 0.000 r(n) ]};
       
    % Relaxation rates, distance and orientation 
    % dependence provided using a function handle
    inter.relaxation={'t1_t2'};
    r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                                   2.00230,1e-3,52,r(n),bet);
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
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=32;                    % Number of XiX DNP blocks (power of 2)
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.shot_spacing=204e-6;
    parameters.addshift=-13e6;
    parameters.el_offs=offsets;

    % Over B1 fields
    for k=1:numel(b1)

        % Set electron nutation frequency
        parameters.irr_powers=b1(k);

        % Run the steady state simulation
        dnp(:,n,k)=powder(spin_system,@xixdnp_steady,parameters,'esr');

    end
     
end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 1 numel(wb1)]),3)/sum(wb1);

% Integrate over the distance distribution, r^2 is the radial part of the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(wr,[1 numel(wr)]),2)/sum((r.^2).*wr);

% Plotting 
figure(); plot(parameters.el_offs/1e6,real(dnp)); 
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Microwave resonance offset, MHz'); kgrid; xlim tight;

end

