% Tests tied Cartesian repetition through the full GRAPE objective.
% Syntax:
%
%                  result=test_grape_repeat()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_grape_repeat()

% Construct a full-time one-spin problem with unequal, noncommuting drift slices
result=new_test_result('kernel/grape_repeat','Tied GRAPE repetition',...
                       'tied controls must preserve the full drift sequence and pull derivatives back.');
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
control.operators={S.y};
control.rho_init={S.x};
control.rho_targ={S.z};
control.pwr_levels=1;
control.pulse_dt=[0.02 0.03 0.04 0.01];
control.drifts={{S.z,-2*S.z,S.x,0.4*S.z}};
control.method='newton';
control.max_iter=4;
control.penalties={'NS'};
control.p_weights=0.2;
control.l_bound=-100;
control.u_bound=100;
control.plotting={};
spin_system=optimcon(spin_system,control);

% Compare the tied objective to evaluation on the complete waveform
waveform=[2 4];
[~,fid,grad,hess]=grape_repeat(waveform,spin_system);
[~,expected]=grape_xy([waveform waveform],spin_system);
result=test_close(result,'expanded fidelity and penalty',fid,expected,0,0,...
                  'physical fidelity and penalties are both evaluated over all repeats');

% Move both copies together and check the physical and penalty derivatives
direction=[-0.7 0.5];
delta=eps^(1/3);
plus=waveform+delta*direction;
minus=waveform-delta*direction;
[~,fid_plus,grad_plus]=grape_xy([plus plus],spin_system);
[~,fid_minus,grad_minus]=grape_xy([minus minus],spin_system);
deriv=reshape(sum(grad.*direction,2),1,[]);
result=test_close(result,'tied directional derivative',deriv,(fid_plus-fid_minus)/(2*delta),...
                  1e-8,1e-8,'a tied perturbation changes the controls at every repeated time');

% Check both Hessian axes against an expanded, tied perturbation
grad_diff=(grad_plus-grad_minus)/(2*delta);
expected=reshape(grad_diff(:,1:2,:)+grad_diff(:,3:4,:),2,[]);
observed=zeros(size(expected));
for n=1:size(hess,3)
    observed(:,n)=hess(:,:,n)*direction.';
end
result=test_close(result,'tied Hessian action',observed,expected,1e-8,1e-8,...
                  'both input and output derivative axes are tied across the full time grid.');

% Optimise the actual short guess through the public optimiser
optimised=fmaxnewton(spin_system,@grape_repeat,waveform);
[~,improved]=grape_repeat(optimised,spin_system);
result=test_true(result,'short guess optimisation',isequal(size(optimised),size(waveform))&&...
                 improved(1)>fid(1),'the public optimiser improves fidelity without untying the repeats.');

end


