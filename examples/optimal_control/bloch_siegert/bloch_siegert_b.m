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
control.method='rbfgs';
control.max_iter=100;
control.plotting={};

% Power levels to sweep (rad/s)
zeeman_frq=abs(spin('13C')*sys.magnet);
pwr_list=zeeman_frq*linspace(1e-3,1,20);

% Preallocate results
inf_no_bss=zeros(size(pwr_list));
inf_bss=zeros(size(pwr_list));

% Common initial guess
guess=randn(2,50)/10;

% Loop over power 
for n=1:numel(pwr_list)

    % Set current power
    control.pwr_levels=pwr_list(n);

    % Set pulse slice duration
    dt=(pi/100)/control.pwr_levels;
    control.pulse_dt=dt*ones(1,50);

    % Optimisation without BSS
    control.bsiegert=false();
    spin_system_no_bss=optimcon(spin_system,control);
    pulse_no_bss=fmaxnewton(spin_system_no_bss,@grape_xy,guess);

    % Optimisation with BSS
    control.bsiegert=true();
    spin_system_bss=optimcon(spin_system,control);
    pulse_bs=fmaxnewton(spin_system_bss,@grape_xy,guess);
    
    % Test the pulse without BSS in reality with BSS
    pulse_no_bss=mat2cell(control.pwr_levels*pulse_no_bss,[1 1]);
    [controls_aug,pulse_aug]=bloch_siegert(spin_system_bss,{Lx,Ly},pulse_no_bss);
    rho=shaped_pulse_xy(spin_system_bss,D,controls_aug,pulse_aug,...
                        control.pulse_dt,[Sx,Sy,Sz],'expv-pwc');
    inf_no_bss(n)=3-real(trace([-Sz Sy Sx]'*rho));

    % Test the pulse with BSS in reality with BSS
    pulse_bss=mat2cell(control.pwr_levels*pulse_bs,[1 1]);
    [controls_aug,pulse_aug]=bloch_siegert(spin_system_bss,{Lx,Ly},pulse_bss);
    rho=shaped_pulse_xy(spin_system_bss,D,controls_aug,pulse_aug,...
                        control.pulse_dt,[Sx,Sy,Sz],'expv-pwc');
    inf_bss(n)=3-real(trace([-Sz Sy Sx]'*rho));
   
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

