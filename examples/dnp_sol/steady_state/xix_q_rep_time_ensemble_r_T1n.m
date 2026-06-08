% Simulation of T1n dependent XiX DNP optimisation of 
% repetition time in the steady state with electron-
% proton distance ensemble.
% 
% Calculation time: minutes
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function xix_q_rep_time_ensemble_r_T1n()

% Nuclear relaxation times, seconds
T1n=[50.0 5.0 0.500 0.050 0.005];

% Get the figure started
kfigure(); hold on; kgrid;
kxlabel('XiX repetition time (ms)');
kylabel('$\langle I_Z \rangle _{\infty}$');
xlim([0 1]); ylim([-0.01e-3 0.7e-3]);

% Plot the curves
for n=1:numel(T1n)
    xix_rep_time_ensemble_r(T1n(n));
end

% Add the legend and save the plot
klegend({'$T_{1n}$ = 50 s', '$T_{1n}$ = 5 s',...
         '$T_{1n}$ = 0.5 s','$T_{1n}$ = 0.05 s',...
         '$T_{1n}$ = 0.005 s'},'Location','Best');
savefig(gcf,'xix_q_rep_time_ensemble_r_T1n.fig');

end

% Simulation for a specific T1n
function xix_rep_time_ensemble_r(T1n)

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
    inter.r1_rates={1e3 1/T1n};
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
    parameters.irr_powers=18e6;              % Electron nutation frequency [Hz]
    parameters.grid='rep_2ang_800pts_sph';
    parameters.pulse_dur=48e-9;              % Pulse duration, seconds
    parameters.nloops=36;                    % Number of XiX DNP blocks
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=-39e6;

    % Over repetition times
    for m=1:numel(rep_time)

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

% Update the figure
plot(1e3*rep_time,-real(dnp)); drawnow();

end

