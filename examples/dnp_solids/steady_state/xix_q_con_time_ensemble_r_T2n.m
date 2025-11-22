% Simulation of T1n dependence of XiX DNP contact 
% curves in the steady state with electron-proton
% distance ensemble.
% 
% Calculation time: hours
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function xix_q_con_time_ensemble_r_T2n()

% Nuclear T2 times, seconds
T2n=[2000e-6 200e-6 20e-6 2e-6 0.2e-6];

% Get the figure started
kfigure(); hold on; kgrid;
kxlabel('XiX contact time ($\mu$s)');
kylabel('$\langle I_Z \rangle _{\infty}$');
xlim tight; ylim([0 14e-4]);

% Plot the curves
for n=1:numel(T2n)
    xix_contact_curve_ensemble_r(T2n(n));
end

% Add the legend and save the plot
klegend({'$T_{2n}$ = 2000 $\mu$s', '$T_{2n}$ = 200 $\mu$s',...
         '$T_{2n}$ = 20 $\mu$s','$T_{2n}$ = 2 $\mu$s',...
         '$T_{2n}$ = 0.2 $\mu$s'},'Location','NorthEast');
savefig(gcf,'xix_q_con_time_ensemble_r_T2n.fig');

end

% Simulation for a specific T2n
function xix_contact_curve_ensemble_r(T2n)

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

% XiX loop count
loop_counts=1:64;

% Distance ensemble
[r,wr]=gaussleg(3.5,20,3);      % Angstrom

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
    inter.r1_rates={1e3 r1n_rate};
    inter.r2_rates={200e3 1/T2n};
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
    parameters.grid='rep_2ang_100pts_oct';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=61e6;
   
    % Over loop counts
    parfor m=1:numel(loop_counts)

        % Localise parameters
        localpar=parameters;

        % Set the number of loops
        localpar.nloops=loop_counts(m);

        % Update the shot spacing
        pulses_dur=2*localpar.nloops*localpar.pulse_dur;
        localpar.shot_spacing=153e-6 - pulses_dur;

        % Run the steady state simulation
        dnp(m,n)=powder(spin_system,@xixdnp_steady,localpar,'esr');

    end

end

% Integrate over the distance distribution, r^2 is the radial part of the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(wr,[1 numel(wr)]),2)/sum((r.^2).*wr);

% Plotting 
contact_times=2*parameters.pulse_dur*loop_counts;
plot(1e6*contact_times,real(dnp)); drawnow;

end

