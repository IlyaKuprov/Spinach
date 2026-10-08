% GRAPE population transfer between the two lowest levels of a Duf-
% fing transmon in the rotating frame of the drive, with the two tra-
% jectory cost terms of Chapter 3 of Yunwei Lu's PhD thesis (North-
% western University, 2026): a running penalty on the population of
% the leakage levels above the qubit subspace, and the fidelity ave-
% raged over the pulse nodes rather than taken at the last node. The
% pulse is short enough for the leakage levels to be within the drive
% bandwidth; the penalty keeps them empty at a small cost in fidelity,
% the time-averaged fidelity brings the population to the target ear-
% ly in the pulse, and the two terms may be combined.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function transmon_grape_leakage()

% Magnet field
sys.magnet=0;

% Four-level transmon on resonance with its drive
sys.isotopes={'T4'};
inter.modes.frqs={0};
inter.modes.anharms={-200e6};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotating frame drift Hamiltonian
H=hamiltonian(assume(spin_system,'labframe'));

% Quadrature control operators
Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1);
ops={(Cr+An)/2,1i*(Cr-An)/2};

% Ground and first excited transmon levels
rho_init=state(spin_system,{'BL1'},{1});
rho_targ=state(spin_system,{'BL2'},{1});

% Projector onto the leakage levels
leak_proj=state(spin_system,{'BL3'},{1})+state(spin_system,{'BL4'},{1});

% Power level and time grid, 6 ns in 24 steps
pwr_level=2*pi*200e6; pulse_dt=0.25e-9*ones(1,24);

% Common control parameters
common.isotopes={'T4'};
common.channels=[1;1];
common.drifts={{H}};
common.operators=ops;
common.rho_init={rho_init};
common.rho_targ={rho_targ};
common.pwr_levels=pwr_level;
common.pulse_dt=pulse_dt;
common.penalties={'NS'};
common.p_weights=0.001;
common.method='lbfgs';
common.max_iter=300;
common.plotting={};

% Trajectory cost term variants
variants={'final fidelity','leakage penalty',...
          'time-averaged fidelity','time-averaged fidelity, leakage penalty'};

% Preallocate fidelity and leakage population trajectories
fid_traj=zeros(4,24); leak_traj=zeros(4,24);

% Random initial guess
rng(1); guess=0.3*randn(2,24);

% Loop over the cost term variants
for v=1:4

    % Trajectory cost terms on top of the common parameters
    control=common;
    if ismember(v,[2 4]), control.traj_pen={leak_proj}; end
    if ismember(v,[3 4]), control.fid_type='average'; end

    % Spinach housekeeping
    spin_system=optimcon(spin_system,control);

    % Run the optimisation, get normalised pulse
    pulse=fmaxnewton(spin_system,@grape_xy,guess);

    % Fidelity and leakage population at every node by direct propagation
    rho=rho_init;
    for n=1:24
        slice_ham=H+pwr_level*(pulse(1,n)*ops{1}+pulse(2,n)*ops{2});
        P=propagator(spin_system,slice_ham,pulse_dt(n));
        rho=P*rho*P'; fid_traj(v,n)=real(trace(rho_targ'*rho));
        leak_traj(v,n)=real(trace(leak_proj'*rho));
    end

    % Report the outcome
    disp([pad(variants{v},40) ' final fidelity ' num2str(fid_traj(v,end),'%.4f') ...
          ', first node above 0.9 fidelity ' int2str(find(fid_traj(v,:)>0.9,1)) ...
          ', mean leakage population ' num2str(mean(leak_traj(v,:)),'%.4f') ...
          ', peak leakage population ' num2str(max(leak_traj(v,:)),'%.4f')]);

end

% Validate the optimisations
if any(fid_traj(:,end)<0.99)
    error('GRAPE optimisation did not converge.');
end
if (mean(leak_traj(2,:))>mean(leak_traj(1,:)))||(mean(leak_traj(4,:))>mean(leak_traj(3,:)))
    error('leakage penalty did not reduce the leakage population.');
end
if find(fid_traj(3,:)>0.9,1)>find(fid_traj(1,:)>0.9,1)
    error('time-averaged fidelity did not bring the target forward.');
end

% Plot the fidelity trajectories
time_axis=1e9*cumsum(pulse_dt);
kfigure(); subplot(1,2,1); plot(time_axis,fid_traj'); kgrid;
kxlabel('time, ns'); kylabel('transfer fidelity');
klegend(variants,'Location','southeast');

% Plot the leakage population trajectories
subplot(1,2,2); plot(time_axis,leak_traj'); kgrid;
kxlabel('time, ns'); kylabel('leakage population');
klegend(variants,'Location','northwest');

end

