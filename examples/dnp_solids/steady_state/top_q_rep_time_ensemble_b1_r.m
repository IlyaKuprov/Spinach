% Simulation of TOP DNP repetition time scan in the steady 
% state with distributions in electron-proton distance and
% microwave B1 field.
% 
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function top_q_rep_time_ensemble_b1_r()

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

% Distance and B1 ensemble
[r,wr]=gaussleg(3.5,20,3);      % Angstrom
[b1,wb1]=gaussleg(10e6,20e6,5); % Hz

% Log spacing for rep. time
rep_time=logspace(-5,-3,30);

% Preallocate equilibrium DNP value array
dnp=zeros([numel(rep_time) numel(r) numel(b1)],'like',1i);

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
    parameters.delay_dur=14e-9;              % Delay duration, seconds
    parameters.nloops=300;                   % Number of TOP DNP blocks
    parameters.addshift=-13e6;
    parameters.el_offs=95e6;

    % Over B1 fields
    for k=1:numel(b1)     
        
        % Set electron nutation frequency
        parameters.irr_powers=b1(k);             
    
        % Over repetition times
        parfor m=1:numel(rep_time)

            % Localise parameters
            localpar=parameters;

            % Set the shot spacing
            pulses_time=localpar.nloops*(localpar.pulse_dur+...
                                         localpar.delay_dur);
            localpar.shot_spacing=rep_time(m)-pulses_time;

            % Run the steady state simulation
            dnp(m,n,k)=powder(spin_system,@topdnp_steady,localpar,'esr');

        end
    
    end

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 1 numel(wb1)]),3)/sum(wb1);

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(wr,[1 numel(wr)]),2)/sum((r.^2).*wr);

% Plotting 
figure(); plot(rep_time*1e3,real(dnp));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Repetition time, ms'); 
kgrid; xlim([0 2]); ylim padded;

% Save for later
savefig(gcf,'top_q_rep_time_ensemble_b1_r.fig');

end

