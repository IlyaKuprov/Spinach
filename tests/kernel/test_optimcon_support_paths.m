% Tests small optimal-control support paths. Syntax:
%
%                    result=test_optimcon_support_paths()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test covers penalty functions, trapezium-product derivatives,
% objective-function collection, and small line-search helper paths.
%
% ilya.kuprov@weizmann.ac.il

function result=test_optimcon_support_paths()

% Announce the test target
fprintf('TESTING: Optimal-control support helper paths\n');

% State the support-path target of the test
result=new_test_result('optimcon/support_paths',...
                       'Optimal-control support helper paths',...
                       'penalty(), trapdiff(), objeval(), and line-search helpers must pass small deterministic checks.');

% Make a minimal quiet Spinach object for low-level helper calls
spin_system=local_spin_system();

% Check the no-penalty path
waveform=[-2.0 -0.5 0.25 2.0 0.75 -1.25;...
           1.5 -1.5 0.10 -0.1 0.50  1.25];
[pen_none,grad_none,hess_none]=penalty(waveform,'none',-1,1);
result=test_close(result,'penalty none term',pen_none,0,0,0,...
                  'the none penalty must have zero value');
result=test_close(result,'penalty none gradient',grad_none,zeros(size(waveform)),0,0,...
                  'the none penalty must have zero gradient');
result=test_close(result,'penalty none Hessian',hess_none,zeros(numel(waveform),numel(waveform)),0,0,...
                  'the none penalty must have zero Hessian');

% Check the norm-square penalty against its closed form
[pen_ns,grad_ns,hess_ns]=penalty(waveform,'NS',-1,1);
pen_ns_ref=sum(waveform(:).^2)/size(waveform,2);
grad_ns_ref=2*waveform/size(waveform,2);
hess_ns_ref=2*eye(numel(waveform))/size(waveform,2);
result=test_close(result,'penalty NS term',pen_ns,pen_ns_ref,1e-14,1e-14,...
                  'NS penalty value is the mean squared waveform norm');
result=test_close(result,'penalty NS gradient',grad_ns,grad_ns_ref,1e-14,1e-14,...
                  'NS penalty gradient is twice the waveform divided by the number of time points');
result=test_close(result,'penalty NS Hessian',hess_ns,hess_ns_ref,1e-14,1e-14,...
                  'NS penalty Hessian is a scaled identity matrix');

% Check the spillout penalty against explicit clipping residuals
[pen_sns,grad_sns,hess_sns]=penalty(waveform,'SNS',-1,1);
spill_hi=max(waveform-1,0);
spill_lo=min(waveform+1,0);
pen_sns_ref=sum(spill_hi(:).^2)+sum(spill_lo(:).^2);
pen_sns_ref=pen_sns_ref/size(waveform,2);
grad_sns_ref=2*(spill_hi+spill_lo)/size(waveform,2);
spill_mask=(abs(spill_hi)+abs(spill_lo)>0);
hess_sns_ref=2*diag(spill_mask(:));
hess_sns_ref=hess_sns_ref/size(waveform,2);
result=test_close(result,'penalty SNS term',pen_sns,pen_sns_ref,1e-14,1e-14,...
                  'SNS penalty value is the mean square distance outside the bounds');
result=test_close(result,'penalty SNS gradient',grad_sns,grad_sns_ref,1e-14,1e-14,...
                  'SNS penalty gradient is twice the clipping residual divided by the number of time points');
result=test_close(result,'penalty SNS Hessian',hess_sns,hess_sns_ref,1e-14,1e-14,...
                  'SNS penalty Hessian is diagonal on points outside the bounds');

% Check the derivative norm-square gradient by finite differences
[pen_dns,grad_dns,hess_dns]=penalty(waveform,'DNS',-10,10);
grad_dns_ref=local_finite_grad(@(x)local_penalty_term(x,size(waveform),'DNS',-10,10),...
                               waveform(:));
result=test_true(result,'penalty DNS finite term',isfinite(pen_dns),...
                 'DNS penalty value must be finite for a small waveform');
result=test_close(result,'penalty DNS gradient',grad_dns(:),grad_dns_ref,1e-6,1e-6,...
                  'DNS penalty gradient must match centred finite differences');
result=test_true(result,'penalty DNS Hessian shape',isequal(size(hess_dns),[numel(waveform) numel(waveform)]),...
                 'DNS penalty Hessian must be square over waveform elements');

% Check the amplitude spillout path on Cartesian controls
amp_waveform=[0.5 2.0 -0.25;...
              0.5 0.0  1.50];
[pen_snsa,grad_snsa,hess_snsa]=penalty(amp_waveform,'SNSA',0,1);
result=test_true(result,'penalty SNSA active',pen_snsa>0,...
                 'SNSA penalty must detect amplitudes above the ceiling');
result=test_true(result,'penalty SNSA gradient shape',isequal(size(grad_snsa),size(amp_waveform)),...
                 'SNSA penalty gradient must have waveform dimensions');
result=test_true(result,'penalty SNSA Hessian shape',isequal(size(hess_snsa),[numel(amp_waveform) numel(amp_waveform)]),...
                 'SNSA penalty Hessian must be square over waveform elements');

% Check trapezium derivative matrices against centred finite differences
S=pauli(2);
Hd={2*pi*3*S.z+0.1*S.x, -2*pi*2*S.z+0.2*S.y};
Hc=S.x+0.3*S.z;
dt=1e-3;
cL=2.0;
cR=-1.5;
[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR);
H_dir_L=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hc*Hd{2}-Hd{2}*Hc);
H_dir_R=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hd{1}*Hc-Hc*Hd{1});
H=(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R;
dc=1e-6;
DL_ref=(expm(-1i*dt*(H+dc*H_dir_L))-expm(-1i*dt*(H-dc*H_dir_L)))/(2*dc);
DR_ref=(expm(-1i*dt*(H+dc*H_dir_R))-expm(-1i*dt*(H-dc*H_dir_R)))/(2*dc);
result=test_close(result,'trapdiff left derivative',DL,DL_ref,1e-10,1e-10,...
                  'left-edge trapezium derivative must match finite differences');
result=test_close(result,'trapdiff right derivative',DR,DR_ref,1e-10,1e-10,...
                  'right-edge trapezium derivative must match finite differences');

% Check objective evaluation for fidelity, gradient, and Hessian collection
data=local_data([2 1]);
x=[0.2; -0.4];
[data_fx,fx]=objeval(x,@local_objective,data,spin_system);
result=test_close(result,'objeval value',fx,local_total_objective(x),1e-14,1e-14,...
                  'objeval must combine fidelity and penalty terms by subtraction');
result=test_true(result,'objeval value count',data_fx.count.fx==1,...
                 'objective-value calls must increment the value counter');
[data_grad,fx_grad,grad]=objeval(x,@local_objective,data,spin_system);
result=test_close(result,'objeval gradient value',fx_grad,local_total_objective(x),1e-14,1e-14,...
                  'gradient calls must return the same objective value');
result=test_close(result,'objeval gradient',grad,local_total_gradient(x),1e-14,1e-14,...
                  'objeval must combine gradient channels by subtraction');
result=test_true(result,'objeval gradient counts',(data_grad.count.fx==1)&&(data_grad.count.gfx==1),...
                 'gradient calls must increment value and gradient counters');
[data_hess,fx_hess,grad_hess,hess]=objeval(x,@local_objective,data,spin_system);
result=test_close(result,'objeval Hessian value',fx_hess,local_total_objective(x),1e-14,1e-14,...
                  'Hessian calls must return the same objective value');
result=test_close(result,'objeval Hessian gradient',grad_hess,local_total_gradient(x),1e-14,1e-14,...
                  'Hessian calls must return the combined gradient');
result=test_close(result,'objeval Hessian',hess,local_total_hessian(),1e-14,1e-14,...
                  'objeval must combine Hessian channels by subtraction');
result=test_true(result,'objeval Hessian counts',(data_hess.count.fx==1)&&...
                 (data_hess.count.gfx==1)&&(data_hess.count.hfx==1),...
                 'Hessian calls must increment all objective counters');

% Check direct line-search condition helpers
spin_system.control.freeze=[];
spin_system.control.ls_c1=1e-2;
spin_system.control.ls_c2=0.9;
spin_system.control.ls_tau1=3;
spin_system.control.ls_tau2=0.1;
spin_system.control.ls_tau3=0.5;
dir=[2; -4];
gfx_0=dir;
gfx_1=[1.6; -3.2];
result=test_true(result,'alpha_conds monotonic',alpha_conds(0,[],1,2,[],[],[],spin_system),...
                 'monotonic line-search test must accept a larger objective value');
result=test_true(result,'alpha_conds Armijo',alpha_conds(1,0.1,-5,-3.2,gfx_0,[],dir,spin_system),...
                 'Armijo test must accept sufficient ascent');
result=test_true(result,'alpha_conds curvature',alpha_conds(2,[],[],[],gfx_0,gfx_1,dir,spin_system),...
                 'strong Wolfe curvature test must accept the reduced directional derivative');
[alpha_cubic,fx_cubic]=cubic_interp(0,1,0,1,0,1,0,-1);
result=test_close(result,'cubic_interp maximiser',alpha_cubic,0.5,1e-14,1e-14,...
                  'Hermite cubic with opposite endpoint slopes has its maximum halfway');
result=test_close(result,'cubic_interp value',fx_cubic,0.25,1e-14,1e-14,...
                  'the corresponding cubic maximum value is one quarter');

% Check bracketing accepts a short ascent step on a concave quadratic
data_line=local_data([2 1]);
x_0=[0; 0];
[data_line,fx_0,gfx_0]=objeval(x_0,@local_line_objective,data_line,spin_system);
dir=gfx_0;
[br_a,br_b,alpha_br,fx_br,gfx_br,next_act,data_line]=bracketing(@local_line_objective,0.1,dir,x_0,...
                                                                 fx_0,gfx_0,data_line,spin_system); %#ok<ASGLU>
result=test_true(result,'bracketing immediate acceptance',strcmp(next_act,'none'),...
                 'a conservative ascent step on a concave quadratic should satisfy Wolfe tests');
result=test_close(result,'bracketing step length',alpha_br,0.1,1e-14,1e-14,...
                  'accepted conservative trial step should be returned unchanged');
result=test_true(result,'bracketing objective increase',fx_br>fx_0,...
                 'accepted line-search step must increase the objective');
result=test_true(result,'bracketing gradient shape',isequal(size(gfx_br),size(gfx_0)),...
                 'bracketing must return a gradient vector with the original dimensions');

% Check sectioning recovers the bracketed maximum of the same quadratic
a.alpha=0;
a.fx=fx_0;
a.gfx=gfx_0;
b.alpha=1;
[~,b.fx,b.gfx]=objeval(x_0+b.alpha*dir,@local_line_objective,local_data([2 1]),spin_system);
[alpha_sec,fx_sec,gfx_sec,exitflag]=sectioning(@local_line_objective,a,b,x_0,fx_0,gfx_0,...
                                               dir,local_data([2 1]),spin_system);
result=test_close(result,'sectioning step length',alpha_sec,0.5,1e-12,1e-12,...
                  'sectioning must locate the one-dimensional quadratic maximum');
result=test_close(result,'sectioning objective',fx_sec,0,1e-12,1e-12,...
                  'the quadratic objective reaches zero at its target');
result=test_close(result,'sectioning gradient',gfx_sec,zeros(size(gfx_sec)),1e-12,1e-12,...
                  'the gradient must vanish at the quadratic maximum');
result=test_true(result,'sectioning exit flag',exitflag==0,...
                 'sectioning should report successful Wolfe acceptance');

end


function spin_system=local_spin_system()

% Build a minimal quiet Spinach object for matrix-level helpers
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.small_matrix=64;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism='zeeman-hilb';

end


function data=local_data(x_shape)

% Build the optimisation workspace fields used by objeval()
data.x_shape=x_shape;
data.count.fx=0;
data.count.gfx=0;
data.count.hfx=0;

end


function term=local_penalty_term(x,wf_size,type,fb,cb)

% Evaluate only the scalar penalty value for finite differences
term=penalty(reshape(x,wf_size),type,fb,cb);

end


function grad=local_finite_grad(fun,x)

% Compute a centred finite-difference gradient
step_size=1e-6;
grad=zeros(size(x));
for n=1:numel(x)
    x_plus=x;
    x_minus=x;
    x_plus(n)=x_plus(n)+step_size;
    x_minus(n)=x_minus(n)-step_size;
    grad(n)=(fun(x_plus)-fun(x_minus))/(2*step_size);
end

end


function [traj_data,fidelity,grad,hess]=local_objective(x,~)

% Define a two-channel objective with a penalty-like second component
traj_data.marker='local_objective';
target=[1; -2];
residual=x-target;
fidelity=zeros(2,1);
fidelity(1)=-sum(residual.^2);
fidelity(2)=0.1*sum(x.^2);
if nargout>2
    grad=zeros([size(x) 2]);
    grad(:,:,1)=-2*residual;
    grad(:,:,2)=0.2*x;
end
if nargout>3
    hess=zeros(numel(x),numel(x),2);
    hess(:,:,1)=-2*eye(numel(x));
    hess(:,:,2)=0.2*eye(numel(x));
end

end


function value=local_total_objective(x)

% Return the scalar objective assembled by objeval()
target=[1; -2];
residual=x-target;
value=-sum(residual.^2)-0.1*sum(x.^2);

end


function grad=local_total_gradient(x)

% Return the vector gradient assembled by objeval()
target=[1; -2];
residual=x-target;
grad=-2*residual-0.2*x;

end


function hess=local_total_hessian()

% Return the Hessian assembled by objeval()
hess=-2.2*eye(2);

end


function [traj_data,fidelity,grad,hess]=local_line_objective(x,~)

% Define a concave quadratic objective for line-search tests
traj_data.marker='local_line_objective';
target=[1; -2];
residual=x-target;
fidelity=-sum(residual.^2);
if nargout>2
    grad=-2*residual;
end
if nargout>3
    hess=-2*eye(numel(x));
end

end


