% Simulation of XiX DNP loop count dependence in the 
% steady state with electron-proton distance ensembles.
%
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% i.kuprov@soton.ac.uk
% guinevere.mathies@uni-konstanz.de

function xix_loop_count_ensemble_r()

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

% Distance ensemble
[r,wr]=gaussleg(3.5,20,3);      % Angstrom

% Number of XiX loops
loop_counts=[1 2 4 8 16 32 64];

% Preallocate equilibrium DNP value array
dnp=zeros([numel(loop_counts) numel(r)],'like',1i);

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
    parameters.irr_powers=17.8e6;            % Electron nutation frequency [Hz]
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.shot_spacing=153e-6;
    parameters.addshift=-13e6;
    parameters.el_offs=61e6;

    % Over loop counts
    for m=1:numel(loop_counts)

        % Set the number of loops
        parameters.nloops=loop_counts(m);

        % Run the steady state simulation
        dnp(m,n)=powder(spin_system,@xixdnp_steady,parameters,'esr');

    end
        
end

% Integrate over the distance distribution, r^2 is the radial part of the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(wr,[1 numel(wr)]),2)/sum((r.^2).*wr);

% Plotting 
contact_times=parameters.pulse_dur*2*loop_counts;
figure(); plot(contact_times*1e6,real(dnp),'-o');
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Contact time, $\mu$s');
kgrid; xlim tight;

end

