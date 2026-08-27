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
control.method='rbfgs';
control.max_iter=100;
control.plotting={};

% Power levels to sweep (rad/s)
zeeman_frq=abs(spin('1H')*sys.magnet);
pwr_list=zeeman_frq*linspace(1e-3,1,20);

% Preallocate results
inf_no_bss=zeros(size(pwr_list));
inf_bss=zeros(size(pwr_list));

% Loop over control power 
for n=1:numel(pwr_list)

    % Set current power
    control.pwr_levels=pwr_list(n);

    % Set pulse slice duration
    dt=(pi/100)/control.pwr_levels;
    control.pulse_dt=dt*ones(1,50);

    % Optimisation without BSS
    control.bsiegert=false();
    spin_system_no_bss=optimcon(spin_system,control);
    pulse_no_bss=fmaxnewton(spin_system_no_bss,@grape_xy,randn(2,50)/10);

    % Optimisation with BSS
    control.bsiegert=true();
    spin_system_bss=optimcon(spin_system,control);
    pulse_bs=fmaxnewton(spin_system_bss,@grape_xy,randn(2,50)/10);

    % Test the pulse without BSS in reality with BSS
    pulse_no_bss=mat2cell(control.pwr_levels*pulse_no_bss,[1 1]);
    [controls_aug,pulse_aug]=bloch_siegert(spin_system_bss,{Lx,Ly},pulse_no_bss);
    rho=shaped_pulse_xy(spin_system_bss,D,controls_aug,pulse_aug,...
                        control.pulse_dt,rho_init,'expv-pwc');
    inf_no_bss(n)=1-real(rho_targ'*rho);

    % Test the pulse with BSS in reality with BSS
    pulse_bss=mat2cell(control.pwr_levels*pulse_bs,[1 1]);
    [controls_aug,pulse_aug]=bloch_siegert(spin_system_bss,{Lx,Ly},pulse_bss);
    rho=shaped_pulse_xy(spin_system_bss,D,controls_aug,pulse_aug,...
                        control.pulse_dt,rho_init,'expv-pwc');
    inf_bss(n)=1-real(rho_targ'*rho);

end

% Compute relative control powers
relative_power=pwr_list./zeeman_frq;

% Plotting
kfigure(); scale_figure([1.0 0.75]); 
plot(relative_power,inf_no_bss); hold on; kgrid;
plot(relative_power,inf_bss); set(gca,'Yscale','log');
kxlabel('relative control power $|\omega_1 / \omega_0|$');
kylabel('terminal infidelity'); ylim padded;
klegend({'GRAPE','GRAPE + BSS'},'Location','East');

end

