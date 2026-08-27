% Bloch-Siegert shift compensation functionality demo. The
% script optimises a universal rotation pulse for a range 
% of resonance offsets. As the control power is increased,
% Bloch-Siegert shift starts to reduce the fidelity unless
% it is correctly accounted for.
%
% Calculation time: minutes.
%
% aditya.dev@weizmann.ac.il

function bloch_siegert_b()

% Magnet field
sys.magnet=28.18;

% 100 non-interacting spins at equal intervals 
% within [-100,+100] ppm chemical shift range 
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
Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2);
Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2);
Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Get the drift Hamiltonian
D=hamiltonian(assume(spin_system,'nmr'));

% Optimal control settings
control.isotopes={'13C'};
control.channels=[1; 1];
control.drifts={{D}};
control.operators={Lx,Ly};
control.rho_init={ Sx Sy Sz};
control.rho_targ={-Sz Sy Sx};
control.method='lbfgs';
control.max_iter=500;
control.tol_x=1e-4;
control.plotting={};

% Power levels to sweep (rad/s)
zeeman_frq=abs(spin('13C')*sys.magnet);
pwr_list=zeeman_frq*linspace(1e-3,1.0,20);

% Preallocate results
fid_a=zeros(size(pwr_list));
fid_b=zeros(size(pwr_list));

% Common initial guess
guess=randn(2,50)/10;

% Loop over control power 
for n=1:numel(pwr_list)

    % Set current power
    control.pwr_levels=pwr_list(n);

    % Set pulse slice duration
    dt=(pi/100)/control.pwr_levels;
    control.pulse_dt=dt*ones(1,50);

    % BSS settings: off, on
    control.bsiegert=false();
    setting_a=optimcon(spin_system,control);
    control.bsiegert=true();  
    setting_b=optimcon(spin_system,control);

    % Optimisation with and without BSS 
    pulse_a=fmaxnewton(setting_a,@grape_xy,guess);
    pulse_b=fmaxnewton(setting_b,@grape_xy,guess);

    % Evaluation for a system with BSS
    [~,fid_a(n)]=ensemble(pulse_a,setting_b);
    [~,fid_b(n)]=ensemble(pulse_b,setting_b);

end

% Compute relative control powers
relative_power=pwr_list./zeeman_frq;

% Plotting
kfigure(); scale_figure([1.0 0.75]); 
plot(relative_power,1-fid_a); hold on; kgrid;
plot(relative_power,1-fid_b); set(gca,'Yscale','log');
kxlabel('relative control power $|\omega_1 / \omega_0|$');
kylabel('terminal infidelity'); ylim padded;
klegend({'GRAPE','GRAPE + BSS'},'Location','East');

end

