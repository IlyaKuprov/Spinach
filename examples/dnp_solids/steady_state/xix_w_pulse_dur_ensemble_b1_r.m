% 2D parameter scan of XiX DNP in the steady state with 
% electron-proton distance and electron Rabi frequency 
% ensembles.
% 
% Calculation time: hours.
% 
% shebha-anandhi.jegadeesan@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function xix_w_pulse_dur_ensemble_b1_r()

% W-band magnet
sys.magnet=3.4;

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
[r,w]=gaussleg(3.5,20,3);
[b1,wb1]=gaussleg(10e6,20e6,5); % Hz

% Electron pulse duration grid, s
pulse_durs=linspace(2e-9,21e-9,200);

% Microwave resonance offsets, Hz
offsets=linspace(-230e6,205e6,101);

% Preallocate steady state DNP array
dnp=zeros([numel(offsets) numel(pulse_durs) ...
           numel(r) numel(b1)],'like',1i);

% Over distances
for n=1:numel(r)

    % Cartesian coordinates
    inter.coordinates={[0.000 0.000 0.000];
                       [0.000 0.000 r(n) ]};
      
    % Relaxation rates, distance and ori. dep. R1n
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
    parameters.coil=state(spin_system,'Lz',2);

    % Experiment parameters
    parameters.spins={'E','1H'};
    parameters.grid='rep_2ang_800pts_sph';
    parameters.phase=pi;                   % Second pulse inverted phase
    parameters.addshift=-33e6;
    parameters.el_offs=offsets;

    % Over B1 fields
    for k=1:numel(b1)

        % Set electron nutation frequency
        parameters.irr_powers=b1(k);

        % Over pulse durations
        parfor m=1:numel(pulse_durs)
    
            % Localise
            localpar=parameters;

            % Set pulse duration
            localpar.pulse_dur=pulse_durs(m);

            % Loop count: contact time / pulse duration per block
    		localpar.nloops=round(360e-9/(2*pulse_durs(m)));

            % Shot spacing: repetition time - contact time
    		localpar.shot_spacing=167e-6 - 360e-9;
    
            % Run the steady state simulation
            dnp(:,m,n,k)=powder(spin_system,@xixdnp_steady,localpar,'esr');

        end
    
    end

end

% Integrate over the B1 field distribution
dnp=sum(dnp.*reshape(wb1,[1 1 1 numel(wb1)]),4)/sum(wb1);

% Integrate over the distance distribution, r^2 is the Jacobian
dnp=sum(dnp.*reshape(r.^2,[1 1 numel(r)]).*reshape(w,[1 1 numel(w)]),3)/sum((r.^2).*w);

% Do the plotting
figure(); imagesc(parameters.el_offs/1e6,pulse_durs*1e9,real(dnp'));
set(gca,'YDir','normal'); kylabel('Pulse duration, ns');
kxlabel('Microwave resonance offset, MHz'); colormap turbo;
kcolourbar('$I_\textrm{z}$ expectation value on $^{1}$H');

% Save for later
savefig(gcf,'xix_w_pulse_dur_ensemble_b1_r.fig');
save('xix_w_pulse_dur_ensemble_b1_r.mat');

end

