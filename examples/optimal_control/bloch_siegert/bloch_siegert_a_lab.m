% Bloch-Siegert shift compensation functionality demo. The
% script optimises a 90-degree pulse (Lz -> Lx) in the ro-
% tating frame, then applies the pulse in the laboratory
% frame with carrier-explicit RF dynamics.
%
% Calculation time: minutes.
%
% aditya.dev@weizmann.ac.il

function bloch_siegert_a_lab()

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

% Run Spinach housekeeping for lab-frame validation
lab=labframe_setup(sys.magnet,'1H');

% Build and normalise the initial state
rho_init=state(spin_system,{'Lz'},{1});
rho_init=rho_init/norm(full(rho_init),2);

% Build and normalise the target state
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
control.max_iter=20;
control.plotting={};

% Power levels to sweep
zeeman_frq=abs(spin('1H')*sys.magnet);
power_grid=1-(1-1e-4)*(1-linspace(0,1,50)).^3;
pwr_list=zeeman_frq*power_grid;

% Preallocate results
inf_no_bss=zeros(size(pwr_list));
inf_bss=zeros(size(pwr_list));

% Use the same initial guess for both optimisations
n_times=600;

% Loop over control power
for n=1:numel(pwr_list)
    guess=rand(2,n_times)/10;

    % Set current power
    control.pwr_levels=pwr_list(n);

    % Set pulse slice duration
    dt=(pi/n_times)/control.pwr_levels;
    control.pulse_dt=dt*ones(1,n_times);

    % Optimise without BSS in the rotating frame
    control.bsiegert=false();
    spin_system_no_bss=optimcon(spin_system,control);
    pulse_no_bss=fmaxnewton(spin_system_no_bss,@grape_xy,guess);

    % Optimise with BSS in the rotating frame
    control.bsiegert=true();
    spin_system_bss=optimcon(spin_system,control);
    pulse_bss=fmaxnewton(spin_system_bss,@grape_xy,guess);

    % Apply the non-BSS pulse in the laboratory frame
    pulse_no_bss=mat2cell(control.pwr_levels*pulse_no_bss,[1 1]);
    inf_no_bss(n)=labframe_inf(lab,pulse_no_bss,control.pulse_dt);

    % Apply the BSS pulse in the laboratory frame
    pulse_bss=mat2cell(control.pwr_levels*pulse_bss,[1 1]);
    inf_bss(n)=labframe_inf(lab,pulse_bss,control.pulse_dt);

end

% Compute relative control powers
relative_power=pwr_list./zeeman_frq;

% Plotting
kfigure();
plot(relative_power,inf_no_bss); hold on; kgrid;
plot(relative_power,inf_bss); set(gca,'Yscale','log');
kxlabel('relative control power $|\omega_1 / \omega_0|$');
kylabel('lab-frame terminal infidelity');
klegend({'GRAPE','GRAPE + BSS'},'Location','East');
ylim padded;

end

function lab=labframe_setup(magnet,isotope)

% Build a one-spin spherical-tensor Spinach system for validation
sys.magnet=magnet;
sys.isotopes={isotope};
sys.output='hush';
sys.parallel={'processes',1};
inter.zeeman.scalar={0.0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get spherical-tensor spin operators and the carrier frequency
lab.Lx=full(operator(spin_system,'Lx',isotope));
lab.Lz=full(operator(spin_system,'Lz',isotope));
lab.Ix=full(state(spin_system,'Lx',isotope));
lab.Ix=lab.Ix/norm(lab.Ix,2);
lab.Iz=full(state(spin_system,'Lz',isotope));
lab.Iz=lab.Iz/norm(lab.Iz,2);
lab.carrier=-spin(isotope)*magnet;
lab.spin_system=spin_system;

end

function infidelity=labframe_inf(lab,pulse,pulse_dt)

% Start from longitudinal magnetisation
rho=lab.Iz;
time_now=0;

% Apply the rotating-frame pulse as a physical lab-frame waveform
for n=1:numel(pulse_dt)

    % Get current Cartesian pulse coefficients
    ux=pulse{1}(n);
    uy=pulse{2}(n);

    % Apply one pulse slice in the lab frame
    [rho,time_now]=labframe_slice(lab,ux,uy,pulse_dt(n),rho,time_now);

end

% Build the rotating-frame target in the laboratory frame
U=expm(-1i*lab.carrier*lab.Lz*time_now);
rho_targ=U*lab.Ix;

% Compute terminal infidelity after lab-frame observation
infidelity=1-real(trace(rho_targ'*rho));

end

function [rho,time_now]=labframe_slice(lab,ux,uy,slice_dt,rho,...
                                       time_now)

% Set carrier-cycle quadrature
nquad=16;
period=2*pi/abs(lab.carrier);
nperiods=floor(slice_dt/period);
rem_dt=slice_dt-nperiods*period;

% Apply complete carrier periods
if nperiods>0
    P=labframe_prop(lab,ux,uy,time_now,period,nquad);
    P=P^nperiods;
    rho=P*rho;
    time_now=time_now+nperiods*period;
end

% Apply the remaining partial carrier period
if rem_dt>0
    nsteps=max(1,ceil(nquad*rem_dt/period));
    P=labframe_prop(lab,ux,uy,time_now,rem_dt,nsteps);
    rho=P*rho;
    time_now=time_now+rem_dt;
end

end

function P=labframe_prop(lab,ux,uy,time_now,duration,nsteps)

% Start the pulse propagator
P=eye(size(lab.Lz));
dt=duration/nsteps;

% Integrate the carrier-explicit RF Hamiltonian using RKMK-DP8
for n=1:nsteps
    L_fun=@(t,rho) lab.carrier*lab.Lz+...
                   2*(ux*cos(lab.carrier*t)-uy*sin(lab.carrier*t))*lab.Lx;
    P=step(lab.spin_system,{L_fun,time_now,'RKMK-DP8'},P,dt);
    time_now=time_now+dt;
end

end


