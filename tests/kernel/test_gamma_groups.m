% Test: kernel/gamma_groups | Gamma-coherent GRAPE
% Tests coherent gamma-group reduction through the GRAPE objective.
% Syntax:
%
%                  result=test_gamma_groups()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_gamma_groups()

% Construct a full-time one-spin problem with opposite drift phases
result=new_test_result('kernel/gamma_groups','Gamma-coherent GRAPE',...
                       'gamma overlaps interfere only inside a physical group.');
sys.isotopes={'1H'};
sys.magnet=0;
sys.output='hush';
sys.parallel={'processes',1};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=basis(create(sys,inter),bas);
S=pauli(2);
control.isotopes={'1H'};
control.channels=1;
control.operators={S.z};
control.rho_init={S.x};
control.rho_targ={S.y};
control.pwr_levels=1;
control.pulse_dt=1;
control.drift_weights=[1 2];
control.drifts={{(pi/2)*S.z},{-(pi/2)*S.z}};
control.fidelity='square';
control.method='lbfgs';
control.penalties={'none'};
control.p_weights=0.2;
control.l_bound=-100;
control.u_bound=100;
control.plotting={};


% Compare opposite phases assigned to one group and to separate groups
control=rmfield(control,'drift_weights');
control.gamma_groups=[3 3];
together=optimcon(spin_system,control);
control.gamma_groups=[3 7];
separate=optimcon(spin_system,control);
[~,cancelled]=grape_xy(0,together);
[~,preserved]=grape_xy(0,separate);
result=test_close(result,'same group cancellation',cancelled(1),0,1e-12,0,...
                  'opposite terminal overlaps cancel before taking the square');
result=test_true(result,'different groups preserve intensity',preserved(1)>0.01,...
                 'separate crystallites add their intensities rather than their complex amplitudes');

% Compare the coherent derivative to a centred difference away from cancellation
control.drift_weights=[1 2];
control.gamma_groups=[3 3];
together=optimcon(spin_system,control);
amp=0.4;
delta=eps^(1/3);
[~,~,grad]=grape_xy(amp,together);
[~,plus]=grape_xy(amp+delta,together);
[~,minus]=grape_xy(amp-delta,together);
result=test_close(result,'coherent objective derivative',grad(1),(plus(1)-minus(1))/(2*delta),...
                  1e-8,1e-8,'the squared mean differentiates the raw complex overlap');

% Keep independent state-target pairs separate even with identical drift labels
control=rmfield(control,'drift_weights');
control.drifts={{(pi/2)*S.z}};
control.gamma_groups=1;
control.rho_init={S.x,S.x};
control.rho_targ={S.y,-S.y};
independent=optimcon(spin_system,control);
[~,separate]=grape_xy(0,independent);
result=test_close(result,'independent state pairs',separate(1),preserved(1),1e-12,0,...
                  'opposite targets in different catalog rows do not interfere');

% Refuse trajectory penalties before building grouped worker data
control.traj_pen={S.z};
try
    optimcon(spin_system,control);
    refused=false;
catch err
    refused=strcmp(err.identifier,'Spinach:optimcon:gammaTrajectoryPenalty');
end
result=test_true(result,'trajectory penalty rejected',refused,...
                 'grouped overlap must not square a trajectory penalty');

% Run the plotting path with a purely imaginary case overlap
control=rmfield(control,'traj_pen');
control.rho_init={S.x};
control.rho_targ={1i*S.y};
control.plotting={'xy_controls'};
cleaner=test_graphics(); %#ok<NASGU>
plotted=optimcon(spin_system,control);
[~,plot_score]=grape_xy(0,plotted);
result=test_close(result,'complex overlap diagnostics',plot_score(1),preserved(1),1e-12,0,...
                  'the plotted scores are real but the complex group overlap is retained');
end


