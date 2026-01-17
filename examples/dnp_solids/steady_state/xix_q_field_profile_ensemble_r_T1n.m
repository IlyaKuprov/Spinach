% Simulation of T1n dependence of XiX DNP field 
% profiles in the steady state with electron-
% proton distance ensemble.
% 
% Calculation time: minutes
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il

function xix_q_field_profile_ensemble_r_T1n()

% Nuclear relaxation times, seconds
T1n=[50.0 5.0 0.500 0.050 0.005];

% Get the figure started
kfigure(); hold on; kgrid;
kxlabel('Microwave resonance offset, MHz'); 
kylabel('$\langle I_Z \rangle _{\infty}$');
xlim tight; ylim([-1.3e-3 1.3e-3]);

% Plot the curves
for n=1:numel(T1n)
    xix_field_profile_ensemble_r(T1n(n));
end

% Add the legend and save the plot
klegend({'$T_{1n}$ = 50 s', '$T_{1n}$ = 5 s',...
         '$T_{1n}$ = 0.5 s','$T_{1n}$ = 0.05 s',...
         '$T_{1n}$ = 0.005 s'},'Location','SouthEast');
savefig(gcf,'xix_q_field_profile_ensemble_r_T1n.fig');

end

% Simulation for a specific T1n
function xix_field_profile_ensemble_r(T1n)

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
    
    % Relaxation rates
    inter.relaxation={'t1_t2'};
    inter.r1_rates={1e3 1/T1n};
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
    parameters.nloops=36;                    % Number of XiX DNP blocks
    parameters.phase=pi;                     % Second pulse inverted phase
    parameters.addshift=-13e6;
    parameters.el_offs=offsets;

    % Calculate shot spacing
    parameters.shot_spacing=204e-6 - 2*parameters.nloops*parameters.pulse_dur;

    % Run the steady state simulation
    dnp(:,n)=powder(spin_system,@xixdnp_steady,parameters,'esr');

end

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(w,[1 numel(w)]),2)/sum((r.^2).*w);

% Plotting 
plot(parameters.el_offs/1e6,real(dnp)); drawnow;

end

