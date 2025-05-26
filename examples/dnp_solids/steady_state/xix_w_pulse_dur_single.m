% 2D parameter scan of XiX DNP in the steady state.
% 
% Calculation time: minutes.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function xix_w_pulse_dur_single()

% W-band magnet
sys.magnet=3.4;

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

% Electron pulse duration grid, s
pulse_durs=linspace(2e-9,21e-9,200);

% Relaxation rates, distance and ori. dep. R1n
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1e-3,52,r_en,bet);
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
parameters.irr_powers=20e6;            % Electron nutation frequency [Hz]
parameters.grid='rep_2ang_800pts_sph';
parameters.nloops=32;
parameters.phase=pi;                   % Second pulse inverted phase
parameters.shot_spacing=167e-6;
parameters.addshift=-33e6;
parameters.el_offs=linspace(-230e6,205e6,101);

% Preallocate steady state DNP array
dnp=zeros([numel(parameters.el_offs) numel(pulse_durs)],'like',1i);

% Get a figure
figure();

% Over pulse durations
for m=1:numel(pulse_durs)

    % Set pulse duration
    parameters.pulse_dur=pulse_durs(m);

    % Run the steady state simulation
    dnp(:,m)=powder(spin_system,@xixdnp_steady,parameters,'esr');

    % Do the plotting
    imagesc(parameters.el_offs/1e6,pulse_durs*1e9,real(dnp'));
    set(gca,'YDir','normal'); kylabel('Pulse duration, ns');
    kcolourbar('$I_\textrm{z}$ expectation value on $^{1}$H');
    kxlabel('Microwave resonance offset, MHz'); 
    colormap turbo; drawnow();

end

% Save for later
savefig(gcf,'xix_w_pulse_dur_single.fig');
save('xix_w_pulse_dur_single.mat');

end

