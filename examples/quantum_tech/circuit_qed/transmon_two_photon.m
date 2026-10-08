% Two-photon transition in a four-level Duffing transmon. The
% drive carrier is placed halfway between the 0-1 and 1-2 tran-
% sition frequencies, where neither single-photon transition is
% resonant, and GRAPE finds a pulse that moves the population
% from the ground state into the second excited state through a
% virtual intermediate state. Model and parameters from Example
% 2 of the GRAPE_SCQ package.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function transmon_two_photon()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T4'};

% Transmon detuning from the carrier and anharmonicity
inter.modes.frqs={100e6};
inter.modes.anharms={-200e6};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Quadrature control operators of the transmon mode
Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1);
Cx=(Cr+An)/2; Cy=1i*(Cr-An)/2;

% Build source and destination states
rho_init=state(spin_system,'BL1',1);
rho_targ=state(spin_system,'BL3',1);
rho_init=rho_init/norm(rho_init,'fro');
rho_targ=rho_targ/norm(rho_targ,'fro');

% Define control parameters
control.isotopes={'T4'};
control.channels=[1;1];
control.drifts={{H}};
control.operators={Cx,Cy};
control.rho_init={rho_init};
control.rho_targ={rho_targ};
control.pwr_levels=2*pi*100e6;
control.pulse_dt=1e-10*ones(1,200);
control.penalties={'NS'};
control.p_weights=0.001;
control.method='lbfgs';
control.max_iter=200;

% Plots during optimisation
control.plotting={'xy_controls'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Smooth deterministic initial guess
guess=[0.1*sin(pi*(1:200)/200); 0.3*cos(pi*(1:200)/200)];

% Run the optimisation, get normalised pulse
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Propagate the optimal pulse and record populations
dts=spin_system.control.pulse_dt; nslices=numel(dts);
rho=rho_init; pops=zeros(4,nslices+1); pops(:,1)=real(diag(rho));
for n=1:nslices
    slice_ham=H+spin_system.control.pwr_levels*(pulse(1,n)*Cx+pulse(2,n)*Cy);
    P=propagator(spin_system,slice_ham,dts(n));
    rho=P*rho*P'; pops(:,n+1)=real(diag(rho));
end

% Validate the two-photon transfer
if pops(3,end)<0.99
    error('two-photon transfer into the second excited state failed.');
end

% Plot the level populations under the optimal pulse
kfigure(); plot(1e9*[0 cumsum(dts)],pops','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('level populations');
ktitle('two-photon transition');
klegend({'$|0\rangle$','$|1\rangle$','$|2\rangle$','$|3\rangle$'},...
        'Location','Best');

end

