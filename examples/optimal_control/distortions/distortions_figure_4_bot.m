% Figure 4 (bottom) from the paper by Rasulov and Kuprov:
%
%            https://arxiv.org/abs/2502.02198
%
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function distortions_figure_4_bot()

% Set the magnetic field
sys.magnet=28.18;

% Put 100 non-interacting spins at equal intervals 
% within the [-100,+100] ppm chemical shift range 
n_spins=100; sys.isotopes=cell(n_spins,1);
for n=1:n_spins
    sys.isotopes{n}='13C';
end
inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins));

% Select a basis set - IK-2 keeps complete basis on each 
% spin in this case, but ignores multi-spin orders
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sx=state(spin_system,'Lx','13C');
Sy=state(spin_system,'Ly','13C');
Sz=state(spin_system,'Lz','13C');
Sx=Sx/norm(full(Sx),2);
Sy=Sy/norm(full(Sy),2);
Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                              % Drift
control.operators={Lx,Ly};                         % Controls
control.rho_init={ Sx Sy Sz};                      % Starting states
control.rho_targ={-Sz Sy Sx};                      % Target states
control.pulse_dt=1e-6*ones(1,125);                 % Pulse interval grid
control.pwr_levels=2*pi*linspace(50e3,70e3,5);     % Power levels, per channel
control.method='lbfgs';                            % Optimisation method
control.penalties={'NS','SNS'};                    % Penalty types
control.p_weights=[0.01 0.10];                     % Penalty weights
control.max_iter=50;                               % Termination condition
control.parallel='ensemble';                       % Parallelisation

% Last 5 slices are dead time
control.freeze=zeros(2,125);
control.freeze(:,121:125)=1;

% Educated guess
load('guess.mat','xy_profile')
guess=xy_profile; guess(:,121:125)=1e-3;

% Amplifier saturation factor ensemble 
sat_factor=linspace(0.9,1.1,5);
sat_factor=sat_factor*mean(control.pwr_levels);
control.distortion=cell(numel(sat_factor),1);
for n=1:numel(sat_factor)
    control.distortion{n}=@(w)amp_tanh(w,sat_factor(n));
end

% Plotting options
control.plotting={'xy_controls','robustness','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run the LBFGS-GRAPE algorithm
xy_profile=fminnewton(spin_system,@grape_xy,guess);

% Benchmark arrays
sat_factor=linspace(0.5,1.5,41); figure(2);
sat_factor=sat_factor*mean(control.pwr_levels);
pwr_levels=2*pi*linspace(40e3,80e3,41);
infidelity=zeros(numel(sat_factor),numel(pwr_levels));

% Saturation factor loop
for n=1:numel(sat_factor)

    % Nutation frequency loop
    parfor k=1:numel(pwr_levels) %#ok<*PFBNS>

        % Apply the power level
        xy_profile_pwr=pwr_levels(k)*xy_profile;
        
        % Apply amplifier saturation distortion
        xy_profile_dist=amp_tanh(xy_profile_pwr,sat_factor(n));

        % Split the channels
        CLx=xy_profile_dist(1,:);
        CLy=xy_profile_dist(2,:);

        % Simulate the pulse
        rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{CLx,CLy},control.pulse_dt,...
                            cell2mat(control.rho_init),'expv-pwc');

        % Get the fidelity
        infidelity(n,k)=1-real(sum(diag(cell2mat(control.rho_targ)'*rho))/3);

    end

    % Plotting
    imagesc(pwr_levels/(2*pi*1e3),sat_factor/mean(control.pwr_levels),log10(infidelity));
    kylabel('amplifer saturation factor');
    kxlabel('RF nutation frequency, kHz');
    kcolourbar('log(infidelity)');
    set(gca,'YDir','normal');
    xline(50,'w--','LineWidth',1);
    xline(70,'w--','LineWidth',1);
    yline(1.1,'w--','LineWidth',1); 
    yline(0.9,'w--','LineWidth',1);
    drawnow;

end

end

