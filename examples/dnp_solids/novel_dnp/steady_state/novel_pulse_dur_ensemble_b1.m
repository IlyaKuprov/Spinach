% 2D parameter scan of NOVEL DNP in the steady state.
% 
% Calculation time: hours.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function novel_pulse_dur_ensemble_b1()

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
sys.enable={'op_cache','ham_cache'};

% B1 ensemble
[b1,wb1]=gaussleg(14e6,16e6,5); % Hz

% Electron pulse duration grid, s
pulse_durs=linspace(2e-9,21e-9,200);

% Microwave resonance offsets, Hz
offsets=linspace(-300e6,300e6,101);

% Relaxation rates, distance and ori. dep. R1n
inter.relaxation={'t1_t2'};
r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature,...
                               2.00230,1e-3,52,r_en,bet);
inter.r1_rates={1000 r1n_rate};
inter.r2_rates={200000 50e3};
inter.rlx_keep='diagonal';
inter.equilibrium='dibari';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Detect the proton
parameters.coil=state(spin_system,'Lz',2);

% Experiment parameters
parameters.spins={'E','1H'};
parameters.grid='rep_2ang_800pts_sph';
parameters.contact_dur=500e-9;           % Pulse duration, seconds
parameters.flippulse=1;                  % 1 for NOVEL, 0 for SE
parameters.flipback=1;                   % 1 for flipback, 0 for no flipback
parameters.addshift=-3.3e6;
parameters.el_offs=offsets;

% Preallocate steady state DNP array
dnp=zeros([numel(offsets) numel(pulse_durs) numel(b1)],'like',1i);

% Over B1 fields
for k=1:numel(b1)

    % Set electron nutation frequency
    parameters.irr_powers=b1(k);

    % Over pulse durations
    for m=1:numel(pulse_durs)

        % Set pulse duration
        parameters.contact_dur=pulse_durs(m);

        % Run the steady state simulation
        dnp(:,m,k)=powder(spin_system,@noveldnp_steady,parameters,'esr');

    end

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 1 numel(wb1)]),3)/sum(wb1);

% Do the plotting
imagesc(parameters.el_offs/1e6,pulse_durs*1e9,real(dnp'));
set(gca,'YDir','normal'); kylabel('Pulse duration, ns');
kxlabel('Microwave resonance offset, MHz');
kcolourbar('$I_\textrm{z}$ expectation value on $^{1}$H');

end

