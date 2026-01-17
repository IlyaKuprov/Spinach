% Simulation of T2e dependence of XiX DNP repetition time 
% profile in the steady state with electron-proton distan-
% ce ensemble.
% 
% Calculation time: minutes
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function xix_q_rep_time_ensemble_r_T2e()

% Electron relaxation times, seconds
T2e=[50e-6 15e-6 5e-6 1.5e-6 0.5e-6];

% Get the figure started
kfigure(); hold on; kgrid;
kxlabel('Repetition time (ms)');
kylabel('$\langle I_Z \rangle _{\infty}$');
xlim tight; ylim([0 2e-3]);

% Plot the curves
for n=1:numel(T2e)
    xix_rep_time_ensemble_r(T2e(n));
end

% Add the legend and save the plot
klegend({'$T_{2e}$ = 50 $\mu$s', '$T_{2e}$ = 15 $\mu$s',...
         '$T_{2e}$ = 5 $\mu$s','$T_{2e}$ = 1.5 $\mu$s',...
         '$T_{2e}$ = 0.5 $\mu$s'},'Location','NorthEast');
savefig(gcf,'xix_q_rep_time_ensemble_r_T2e.fig');

end

% Simulation for a specific T2e
function xix_rep_time_ensemble_r(T2e)

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
sys.disable={'hygiene'};

% Distance ensemble
[r,wr]=gaussleg(3.5,20,3);      % Angstrom

% Log spacing for rep. time
rep_time=logspace(-5,-3,30);

% Preallocate equilibrium DNP value array
dnp=zeros([numel(rep_time) numel(r)],'like',1i);

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
    inter.r2_rates={1/T2e 50e3};
    inter.rlx_keep='diagonal';
    inter.equilibrium='dibari';
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Detect the proton
    parameters.coil=state(spin_system,'Lz','1H');

    % Experiment parameters
    parameters.spins={'E','1H'};
    parameters.irr_powers=18e6;              % Electron nutation frequency [Hz]
    parameters.grid='rep_2ang_800pts_oct';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=36;                    % Number of XiX DNP blocks
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=-39e6;

    % Over repetition times
    parfor m=1:numel(rep_time)

        % Localise parameters
        localpar=parameters;

        % Set the shot spacing
        pulses_time=2*localpar.nloops*localpar.pulse_dur;
        localpar.shot_spacing=rep_time(m)-pulses_time;

        % Run the steady state simulation
        dnp(m,n)=powder(spin_system,@xixdnp_steady,localpar,'esr');

    end
   
end

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(wr,[1 numel(wr)]),2)/sum((r.^2).*wr);

% Plotting 
plot(rep_time*1e3,-real(dnp)); drawnow;

end

