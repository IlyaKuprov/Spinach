% Simulation of NOVEL DNP repetition time scan in the 
% steady state.
% 
% Calculation time: seconds.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function novel_x_rep_time_single()

% X-band magnet
sys.magnet=0.34;

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

% Relaxation rates, distance and orientation
% dependence provided using a function handle
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1e-3,26,r_en,bet);
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
parameters.irr_powers=15e6;                        % Electron nutation frequency [Hz]
parameters.contact_dur=500e-9;                     % Contact pulse duration, seconds
parameters.pulse_dur=1/(4*parameters.irr_powers);  % 90-degree pulse duration, seconds
parameters.flippulse=1;                            % 1 for NOVEL, 0 for Solid Effect
parameters.addshift=-3.3e6;
parameters.el_offs=0e6;

% Log spacing for rep. time
rep_time=logspace(-4,-2,30);

% Preallocate equilibrium DNP value arrays
dnp_noflip=zeros(size(rep_time),'like',1i);
dnp_flip=zeros(size(rep_time),'like',1i);

% Over repetition times
parfor m=1:numel(rep_time)

    % Localise parameters
    localpar=parameters;

    % Set the shot spacing, no flipback pulse
    localpar.flipback=0;
    pulses_time=localpar.pulse_dur+...
                localpar.contact_dur;
    localpar.shot_spacing=rep_time(m)-pulses_time;

    % Run the steady state simulation
    dnp_noflip(m)=powder(spin_system,@noveldnp_steady,localpar,'esr');

    % Set the shot spacing, flipback pulse
    localpar.flipback=1;
    pulses_time=2*localpar.pulse_dur+...
                  localpar.contact_dur;
    localpar.shot_spacing=rep_time(m)-pulses_time;

    % Run the steady state simulation
    dnp_flip(m)=powder(spin_system,@noveldnp_steady,localpar,'esr');

end

% Plotting 
figure(); plot(rep_time*1e3,real(dnp_noflip));
hold on; plot(rep_time*1e3,real(dnp_flip));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');
klegend({'without flipback','with flipback'});
kxlabel('Repetition time, ms'); 
kgrid; xlim([-1 11]); ylim padded;

% Save for later
savefig(gcf,'novel_x_rep_time_single.fig');
save('novel_x_rep_time_single.mat');

end

