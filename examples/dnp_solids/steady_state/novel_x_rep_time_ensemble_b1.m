% Simulation of NOVEL DNP repetition time scan in the steady 
% state with distributions in microwave B1 field.
% 
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function novel_x_rep_time_ensemble_b1()

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

% B1 ensemble
[b1,wb1]=gaussleg(14e6,16e6,5); % Hz

% Log spacing for rep. time
rep_time=logspace(-4,-2,30);

% Preallocate equilibrium DNP value array
dnp_noflip=zeros([numel(rep_time) numel(b1)],'like',1i);
dnp_flip=zeros([numel(rep_time) numel(b1)],'like',1i);

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
parameters.contact_dur=500e-9;                     % Contact pulse duration, seconds
parameters.flippulse=1;                            % 1 for NOVEL, 0 for Solid Effect
parameters.addshift=-3.3e6;
parameters.el_offs=0e6;

% Over B1 fields
for k=1:numel(b1)

    % Set electron nutation frequency
    parameters.irr_powers=b1(k);

    % recalculate 90-degree pulse duration
    parameters.pulse_dur=1/(4*parameters.irr_powers);  

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
        dnp_noflip(m,k)=powder(spin_system,@noveldnp_steady,localpar,'esr');

        % Set the shot spacing, flipback pulse
        localpar.flipback=1;
        pulses_time=2*localpar.pulse_dur+...
                      localpar.contact_dur;
        localpar.shot_spacing=rep_time(m)-pulses_time;

        % Run the steady state simulation
        dnp_flip(m,k)=powder(spin_system,@noveldnp_steady,localpar,'esr');

    end

end

% Integrate over the B1 field distribution
dnp_noflip=sum(dnp_noflip.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);
dnp_flip=sum(dnp_flip.*reshape(wb1,[1 numel(wb1)]),2)/sum(wb1);

% Plotting 
kfigure(); plot(rep_time*1e3,real(dnp_noflip));
hold on; plot(rep_time*1e3,real(dnp_flip));
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');
klegend({'without flipback','with flipback'});
kxlabel('Repetition time, ms'); 
kgrid; xlim([-1 11]); ylim padded;

% Save for later
savefig(gcf,'novel_x_rep_time_ensemble_b1.fig');
save('novel_x_rep_time_ensemble_b1.mat');

end

