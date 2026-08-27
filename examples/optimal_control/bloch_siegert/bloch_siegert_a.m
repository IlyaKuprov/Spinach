% Bloch-Siegert shift compensation functionality demo. The
% script optimises a 90-degree pulse (Lz -> Lx) for a sing-
% le spin on resonance. As the control power is increased,
% Bloch-Siegert shift starts to reduce the fidelity unless
% it is correctly accounted for.
%
% Calculation time: minutes.
%
% aditya.dev@weizmann.ac.il

function bloch_siegert_a()

% Set magnetic field 
sys.magnet=14.1;

% Set isotopes
sys.isotopes={'1H'};

% Set interactions 
inter.zeeman.scalar={0.0};

% Set basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Build and normalise the initial state (Lz)
rho_init=state(spin_system,{'Lz'},{1});
rho_init=rho_init/norm(full(rho_init),2);

% Build and normalise the target state (Lx)
rho_targ=state(spin_system,{'Lx'},{1});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get the control operators
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');

% Build the drift Hamiltonian
D=hamiltonian(assume(spin_system,'nmr'));

% Optimal control settings
control.isotopes={'1H'};
control.channels=[1; 1];
control.drifts={{D}};
control.operators={Lx,Ly};
control.rho_init={rho_init};
control.rho_targ={rho_targ};
control.method='lbfgs';
control.max_iter=500;
control.tol_x=1e-4;
control.plotting={};

% Power levels to sweep (rad/s)
zeeman_frq=abs(spin('1H')*sys.magnet);
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

