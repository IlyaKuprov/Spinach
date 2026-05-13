% Tests trajectory-analysis dynamic front-end kernels. Syntax:
%
%                    result=test_dynamic_trajectory_frontends()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises trajan() plotting branches and trajsimil() scoring
% branches on a compact two-spin spherical-tensor trajectory.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_trajectory_frontends()

% Announce the test target
fprintf('TESTING: Dynamic trajectory analysis front ends\n');

% State the dynamic trajectory target of the test
result=new_test_result('kernel/dynamic_trajectory_frontends',...
                       'Dynamic trajectory analysis front ends',...
                       'trajectory plotting and similarity helpers must expose deterministic branch outputs.');

% Force invisible figures during plotting checks
old_visibility=get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
cleaner=onCleanup(@()local_cleanup(old_visibility));

% Build the trajectory used by plotting and similarity checks
[spin_system,traj]=local_test_trajectory();

% Check all trajan() property branches
result=local_test_trajan(result,spin_system,traj);

% Check all trajsimil() scoring families
result=local_test_trajsimil(result,spin_system,traj);

end


function result=local_test_trajan(result,spin_system,traj)

% Check correlation-order analysis with an explicit time axis
fig=figure('Visible','off');
time_axis=[0.0 1.0 2.0];
trajan(spin_system,traj,'correlation_order',time_axis);
line_obj=findobj(gca,'Type','line');
result=test_true(result,'trajan correlation-order lines',numel(line_obj)==2,...
                 'trajan() correlation-order mode must plot one-spin and two-spin orders');
result=test_close(result,'trajan explicit time axis',sort(line_obj(1).XData),time_axis,1e-14,1e-14,...
                  'trajan() must use the supplied time axis for plotted trajectory points');
close(fig);

% Check coherence-order analysis
fig=figure('Visible','off');
trajan(spin_system,traj,'coherence_order');
line_obj=findobj(gca,'Type','line');
result=test_true(result,'trajan coherence-order lines',numel(line_obj)==5,...
                 'trajan() coherence-order mode must plot orders from -2 to +2 for two spins');
close(fig);

% Check total population touching each spin
fig=figure('Visible','off');
trajan(spin_system,traj,'total_each_spin');
line_obj=findobj(gca,'Type','line');
result=test_true(result,'trajan total-each-spin lines',numel(line_obj)==2,...
                 'trajan() total_each_spin mode must plot one trace per spin');
close(fig);

% Check local population on each spin
fig=figure('Visible','off');
trajan(spin_system,traj,'local_each_spin');
line_obj=findobj(gca,'Type','line');
result=test_true(result,'trajan local-each-spin lines',numel(line_obj)==2,...
                 'trajan() local_each_spin mode must plot one local trace per spin');
close(fig);

% Check Zeeman level-population branch
fig=figure('Visible','off');
trajan(spin_system,traj,'level_populations');
line_obj=findobj(gca,'Type','line');
result=test_true(result,'trajan level-population lines',numel(line_obj)==4,...
                 'trajan() level_populations mode must plot one trace per Zeeman energy level');
close(fig);

end


function result=local_test_trajsimil(result,spin_system,traj)

% Use identical trajectories for exact similarity references
traj_ref=traj;
score_obs=trajsimil(spin_system,traj,traj_ref,'RSP');
score_ref=sum(conj(traj).*traj_ref,1);
result=test_close(result,'trajsimil running scalar product',score_obs,score_ref,1e-14,1e-14,...
                  'trajsimil() RSP mode must return pointwise trajectory scalar products');

% Check running difference norm for identical trajectories
score_obs=trajsimil(spin_system,traj,traj_ref,'RDN');
result=test_close(result,'trajsimil running difference norm',score_obs,ones(1,size(traj,2)),1e-14,1e-14,...
                  'trajsimil() RDN mode must return unit similarity for identical trajectories');

% Check sign-grouped and broad-state-grouped difference norms
score_obs=trajsimil(spin_system,traj,traj_ref,'SG-RDN');
result=test_close(result,'trajsimil sign-grouped norm',score_obs,ones(1,size(traj,2)),1e-14,1e-14,...
                  'trajsimil() SG-RDN mode must preserve identity under sign grouping');
score_obs=trajsimil(spin_system,traj,traj_ref,'BSG-RDN');
result=test_close(result,'trajsimil broad-grouped norm',score_obs,ones(1,size(traj,2)),1e-14,1e-14,...
                  'trajsimil() BSG-RDN mode must preserve identity under broad state grouping');

end


function [spin_system,traj]=local_test_trajectory()

% Build a two-spin spherical-tensor Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={1.0,2.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10.0;
inter.coupling.scalar{2,2}=0.0;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Assemble a deterministic mixed-order trajectory from physical states
rho_a=unit_state(spin_system)+...
      0.2*state(spin_system,'Lz','1H')+...
      0.1*state(spin_system,'Lz','13C')+...
      state(spin_system,'Lx','1H')+...
      0.5*state(spin_system,'Ly','13C')+...
      0.25*state(spin_system,{'L+','L-'},{1,2});
traj=[rho_a 2*rho_a 0.5*rho_a];

end


function local_cleanup(old_visibility)

% Restore figure state after success or failure
close all force;
set(groot,'defaultFigureVisible',old_visibility);

end


