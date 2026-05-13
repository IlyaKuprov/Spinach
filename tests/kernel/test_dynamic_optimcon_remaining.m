% Tests remaining dynamic optimal-control helper paths. Syntax:
%
%              result=test_dynamic_optimcon_remaining()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test covers the remaining optimcon helpers with small deterministic
% fixtures: waveform distortions, FIR kernel estimation, quasi-Newton
% updates, Hessian handling, waveform utilities, GRAPE wrappers, Liouville
% GRAPE derivatives, TGRAPE duration gradients, fmaxnewton zero-iteration
% handling, and diagnostic plotting smoke paths.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_optimcon_remaining()

% Announce the test target
fprintf('TESTING: Remaining optimal-control dynamic helpers\n');

% State the dynamic optimal-control target of the test
result=new_test_result('kernel/dynamic_optimcon_remaining',...
                       'Remaining optimal-control dynamic helpers',...
                       ['Remaining optimcon helper functions must pass '...
                        'small deterministic regression checks.']);

% Ensure that optimcon and ensemble can use a process-pool ValueStore
local_ensure_pool();

% Run independent groups of small checks
result=local_check_distortions(result);
result=local_check_quasi_newton(result);
result=local_check_wave_utils(result);
result=local_check_grape_family(result);

end


function result=local_check_distortions(result)

% Make a non-degenerate two-channel waveform
waveform=[0.10  0.40 -0.20  1.00;...
          0.20 -0.30  0.50 -0.50];

% Check the identity distortion and its Jacobian
[wave_none,J_none]=no_dist(waveform);
result=test_close(result,'no_dist identity',wave_none,waveform,0,0,...
                  'no_dist must leave the waveform unchanged');
result=test_close(result,'no_dist Jacobian',J_none,speye(numel(waveform)),0,0,...
                  'no_dist must return the vectorised identity Jacobian');

% Check non-orthogonal channel mixing against its closed form
xy_ang=60;
[wave_non,J_non]=non_orth(waveform,xy_ang);
wave_non_ref=[waveform(1,:)+cosd(xy_ang)*waveform(2,:);...
              sind(xy_ang)*waveform(2,:)];
result=test_close(result,'non_orth waveform',wave_non,wave_non_ref,1e-14,1e-14,...
                  'non_orth must tilt the Y channel by the requested angle');
result=local_check_jacobian(result,'non_orth Jacobian',...
                            @(x)non_orth(reshape(x,size(waveform)),xy_ang),...
                            waveform,J_non,1e-9);

% Check FIR filtering against explicit complex convolution
ker=[1.00; 0.25+0.10i; -0.15];
[wave_fir,J_fir]=firf(waveform,ker);
wave_fir_ref=local_firf_ref(waveform,ker);
result=test_close(result,'firf waveform',wave_fir,wave_fir_ref,1e-14,1e-14,...
                  'firf must match truncated causal complex convolution');
result=local_check_jacobian(result,'firf Jacobian',...
                            @(x)firf(reshape(x,size(waveform)),ker),...
                            waveform,J_fir,1e-9);

% Check the single-pole filter recurrence and Jacobian
pole=0.35+0.15i;
[wave_spf,J_spf]=spf(waveform,pole);
wave_spf_ref=local_spf_ref(waveform,pole);
result=test_close(result,'spf waveform',wave_spf,wave_spf_ref,1e-14,1e-14,...
                  'spf must match the causal single-pole recurrence');
result=local_check_jacobian(result,'spf Jacobian',...
                            @(x)spf(reshape(x,size(waveform)),pole),...
                            waveform,J_spf,1e-9);

% Check the single-zero filter recurrence and Jacobian
zero=0.20+0.10i;
[wave_szf,J_szf]=szf(waveform,zero);
wave_szf_ref=local_szf_ref(waveform,zero);
result=test_close(result,'szf waveform',wave_szf,wave_szf_ref,1e-14,1e-14,...
                  'szf must match the causal single-zero recurrence');
result=local_check_jacobian(result,'szf Jacobian',...
                            @(x)szf(reshape(x,size(waveform)),zero),...
                            waveform,J_szf,1e-9);

% Check hyperbolic-tangent amplifier compression
sat_lvl=2.0;
[wave_tanh,J_tanh]=amp_tanh(waveform,sat_lvl);
wave_tanh_ref=local_amp_tanh_ref(waveform,sat_lvl);
amp_tanh_val=sqrt(sum(wave_tanh.^2,1));
result=test_close(result,'amp_tanh waveform',wave_tanh,wave_tanh_ref,1e-14,1e-14,...
                  'amp_tanh must preserve phase and apply tanh amplitude compression');
result=test_true(result,'amp_tanh saturation',all(amp_tanh_val<sat_lvl),...
                 'amp_tanh amplitudes must remain below the saturation level');
result=local_check_jacobian(result,'amp_tanh Jacobian',...
                            @(x)amp_tanh(reshape(x,size(waveform)),sat_lvl),...
                            waveform,J_tanh,1e-8);

% Check root-sigmoidal amplifier compression
shape=4;
[wave_root,J_root]=amp_root(waveform,sat_lvl,shape);
wave_root_ref=local_amp_root_ref(waveform,sat_lvl,shape);
amp_root_val=sqrt(sum(wave_root.^2,1));
result=test_close(result,'amp_root waveform',wave_root,wave_root_ref,1e-14,1e-14,...
                  'amp_root must preserve phase and apply root-sigmoid compression');
result=test_true(result,'amp_root monotonic compression',...
                 all(amp_root_val<=sqrt(sum(waveform.^2,1))),...
                 'amp_root must not increase radial amplitudes');
result=local_check_jacobian(result,'amp_root Jacobian',...
                            @(x)amp_root(reshape(x,size(waveform)),sat_lvl,shape),...
                            waveform,J_root,1e-8);

% Check causal FIR kernel estimation with all solver paths
x=[1; 0; -1; 2; 1; -2; 3];
h_ref=[0.60; -0.20; 0.10];
y=conv(x,h_ref); y=y(1:numel(x));
methods={'backslash','pinv','svd'};
for n=1:numel(methods)
    h_est=kernelest(x,y,numel(h_ref),methods{n},'causal');
    result=test_close(result,['kernelest ' methods{n}],h_est,h_ref,1e-12,1e-12,...
                      'kernelest must recover a full-rank causal FIR kernel');
end
h_tikh=kernelest(x,y,numel(h_ref),'tikh','causal',1e-12);
result=test_close(result,'kernelest tikh',h_tikh,h_ref,1e-8,1e-8,...
                  'Tikhonov-regularised kernel estimation must approach the exact kernel');

% Check same-alignment kernel estimation
conv_mat=toeplitz([x; zeros(numel(h_ref)-1,1)],...
                  [x(1) zeros(1,numel(h_ref)-1)]);
row_start=floor(numel(h_ref)/2)+1;
y_same=conv_mat(row_start:(row_start+numel(x)-1),:)*h_ref;
h_same=kernelest(x,y_same,numel(h_ref),'backslash','same');
result=test_close(result,'kernelest same align',h_same,h_ref,1e-12,1e-12,...
                  'same-aligned kernel estimation must use the central convolution rows');

end


function result=local_check_quasi_newton(result)

% Check a good one-step BFGS update against a known diagonal Hessian
H=bfgs_upd([], [1; 0], [-2; 0]);
result=test_close(result,'bfgs_upd first good pair',H,2*eye(2),1e-14,1e-14,...
                  'a single exact curvature pair must initialise the Hessian scale');

% Check bad-pair safeguards
H_bad=bfgs_upd([], [1; 0], [2; 0]);
result=test_close(result,'bfgs_upd first bad pair',H_bad,eye(2),1e-14,1e-14,...
                  'a first curvature pair with the wrong sign must fall back to identity');
H_unchanged=bfgs_upd(3*eye(2), [1; 0], [2; 0]);
result=test_close(result,'bfgs_upd bad pair unchanged',H_unchanged,3*eye(2),1e-14,1e-14,...
                  'a later bad curvature pair must leave the Hessian unchanged');

% Check full BFGS history reconstruction
x_hist=[1 0; 0 1];
g_hist=[-2 0; 0 -3];
H_hist=bfgs(x_hist,g_hist,[2; 3]);
result=test_close(result,'bfgs diagonal Hessian',H_hist,diag([2 3]),1e-14,1e-14,...
                  'orthogonal curvature pairs must reconstruct a diagonal Hessian');

% Check LBFGS inverse-Hessian action on the same history
lbfgs_dir=lbfgs(x_hist,g_hist,[2; 3]);
result=test_close(result,'lbfgs inverse action',lbfgs_dir,[1; 1],1e-14,1e-14,...
                  'LBFGS must apply the inverse of the reconstructed Hessian');

% Check Hessian reordering against explicit tensor permutation
hess=reshape(1:36,6,6);
hess_ref=reshape(permute(reshape(hess,[2 3 2 3]),[2 1 4 3]),[6 6]);
hess_new=hess_reorder(hess,2,3);
result=test_close(result,'hess_reorder permutation',hess_new,hess_ref,0,0,...
                  'hess_reorder must swap control-first and time-first ordering');

% Check RFO Hessian regularisation on an indefinite Hessian
spin_system.control.reg_alpha=1;
spin_system.control.reg_phi=2;
spin_system.control.reg_max_iter=8;
spin_system.control.reg_max_cond=1e4;
data.count.rfo=0;
[H_reg,data]=hessreg(spin_system,[-1 0; 0 2],[1; 0.5],data);
[~,chol_flag]=chol(H_reg);
result=test_true(result,'hessreg positive definite',chol_flag==0,...
                 'hessreg must shift an indefinite Hessian to positive definiteness');
result=test_true(result,'hessreg counter',data.count.rfo>0,...
                 'hessreg must count the RFO iterations that it performs');
result=test_true(result,'hessreg finite',all(isfinite(H_reg),'all'),...
                 'hessreg must return finite Hessian elements');

end


function result=local_check_wave_utils(result)

% Check frequency-amplitude-phase-time conversion on an explicit grid
fapt={[1.00 2.00 pi/4 0.00 1.00],...
      [0.50 1.00 0.00 0.25 0.75]};
time_grid=0:0.25:1;
[wave,dt,grid_out]=fapt2sfo(fapt,time_grid);
wave_ref=local_fapt_ref(fapt,time_grid);
result=test_close(result,'fapt2sfo waveform',wave,wave_ref,1e-14,1e-14,...
                  'fapt2sfo must add all active rotating components on the grid');
result=test_true(result,'fapt2sfo supplied dt',isempty(dt),...
                 'fapt2sfo must leave dt empty when a grid is supplied');
result=test_close(result,'fapt2sfo grid passthrough',grid_out,time_grid,0,0,...
                  'fapt2sfo must preserve a user-supplied time grid');

% Check instantaneous frequency for a quadratic phase signal
sample_dt=1e-3;
time_axis=0:sample_dt:(6*sample_dt);
base_freq=3;
chirp_rate=4;
signal=exp(1i*2*pi*(base_freq*time_axis+0.5*chirp_rate*time_axis.^2));
freq_ref=base_freq+chirp_rate*time_axis;
freq=inst_freq(signal,sample_dt);
result=test_close(result,'inst_freq chirp',freq,freq_ref,1e-10,1e-10,...
                  'inst_freq must recover the exact frequency of a quadratic phase');

% Check drift extraction through a minimal context callback
spin_system.bas.basis=speye(2);
parameters.marker=true;
[drift_cells,spc_dim]=drifts(spin_system,@local_context,parameters,'nmr');
H=sparse(diag(1:4));
result=test_true(result,'drifts ensemble count',numel(drift_cells)==2,...
                 'drifts must preserve the number of context ensemble members');
result=test_close(result,'drifts first member',drift_cells{1}{1},H+3i*H,1e-14,1e-14,...
                  'drifts must add Hamiltonian, relaxation, and kinetics components');
result=test_close(result,'drifts hydro member',drift_cells{2}{1},2*H+6i*H,1e-14,1e-14,...
                  'drifts must include hydrodynamics when supplied by the context');
result=test_true(result,'drifts spatial dimension',spc_dim==2,...
                 'drifts must infer the classical subspace dimension');

% Check trapezium-product auxiliary matrices
S=pauli(2);
drift_pair={0.10*S.z, -0.20*S.z};
controls={S.x, S.y};
cc_comm=cell(2,2);
cc_comm_idx=false(2,2);
for n=1:2
    for k=1:2
        cc_comm{n,k}=controls{n}*controls{k}-controls{k}*controls{n};
    end
end
time_step=1e-3;
c_left=[0.5; -0.3];
c_right=[0.2; 0.4];
[aux_l,aux_r]=aux_mat(drift_pair,controls,cc_comm_idx,cc_comm,...
                      time_step,c_left,c_right,1);
G=local_trap_generator(drift_pair,controls,time_step,c_left,c_right);
DL=local_trap_dir(drift_pair,controls,cc_comm,cc_comm_idx,time_step,c_left,c_right,1,'left');
DR=local_trap_dir(drift_pair,controls,cc_comm,cc_comm_idx,time_step,c_left,c_right,1,'right');
result=test_close(result,'aux_mat left block',aux_l(1:2,1:2),G,1e-14,1e-14,...
                  'aux_mat must put the trapezium generator on the diagonal');
result=test_close(result,'aux_mat left derivative',aux_l(1:2,3:4),DL,1e-14,1e-14,...
                  'aux_mat left block must match the left directional derivative');
result=test_close(result,'aux_mat right derivative',aux_r(1:2,3:4),DR,1e-14,1e-14,...
                  'aux_mat right block must match the right directional derivative');
[aux_l3,aux_r3]=aux_mat(drift_pair,controls,cc_comm_idx,cc_comm,...
                        time_step,c_left,c_right,1,2);
result=test_true(result,'aux_mat 3x3 sizes',...
                 isequal(size(aux_l3),[6 6])&&isequal(size(aux_r3),[6 6]),...
                 'aux_mat must build 3x3 block matrices for mixed derivatives');

% Check a lightweight diagnostic-plotting branch without showing figures
plot_system.control.plotting={'xy_controls'};
plot_system.control.pulse_dt=[0.01 0.02 0.03];
plot_system.control.pwr_levels=1;
plot_system.control.integrator='rectangle';
plot_system.control.l_bound=-10;
plot_system.control.u_bound=10;
plot_wave=[1 2 3; -1 0 1];
ctrl_trajan(plot_system,plot_wave,{struct()},0.5);
result.messages{end+1}='PASS: ctrl_trajan xy_controls -- plotting smoke path completed offscreen';

end


function result=local_check_grape_family(result)

% Build a tiny Hilbert-space optimal-control problem
spin_system=local_hilb_control_system();
waveform=[4.0  5.0 6.0;...
          1.0 -2.0 3.0];

% Check direct ensemble and Cartesian wrapper consistency
[traj_ens,fid_ens,grad_ens]=ensemble(waveform,spin_system); %#ok<ASGLU>
[traj_xy,fid_xy,grad_xy]=grape_xy(waveform,spin_system); %#ok<ASGLU>
result=test_close(result,'ensemble grape_xy fidelity',fid_xy(1),fid_ens,1e-13,1e-13,...
                  'grape_xy must report the ensemble fidelity in its first channel');
result=test_close(result,'ensemble grape_xy gradient',grad_xy(:,:,1),grad_ens,1e-8,1e-8,...
                  'grape_xy must report the ensemble gradient in its first channel');
result=test_true(result,'grape_xy penalty shape',size(grad_xy,3)==2,...
                 'grape_xy must append one gradient channel per penalty term');

% Check curvilinear wrapper for an identity coordinate transform
[~,fid_curv,grad_curv]=grape_curv(waveform,@local_u2x,@local_dx_du,spin_system);
result=test_close(result,'grape_curv identity fidelity',fid_curv,fid_xy,1e-13,1e-13,...
                  'identity curvilinear coordinates must preserve the fidelity channels');
result=test_close(result,'grape_curv identity gradient',grad_curv,grad_xy,1e-8,1e-8,...
                  'identity curvilinear coordinates must preserve gradient channels');

% Check phase wrapper against Cartesian waveform and finite differences
phi_profile=[0.2 -0.4 0.7];
amp_profile=spin_system.control.amplitudes;
wave_phase=[amp_profile.*cos(phi_profile);...
            amp_profile.*sin(phi_profile)];
[~,fid_phase,grad_phase]=grape_phase(phi_profile,spin_system);
[~,fid_cart]=grape_xy(wave_phase,spin_system);
grad_phase_ref=local_phase_grad(phi_profile,spin_system);
result=test_close(result,'grape_phase fidelity',fid_phase,fid_cart,1e-13,1e-13,...
                  'grape_phase must be equivalent to its Cartesian waveform');
result=test_close(result,'grape_phase gradient',grad_phase(:,:,1),grad_phase_ref,1e-6,1e-6,...
                  'grape_phase gradient must match centred finite differences');

% Check a zero-iteration fmaxnewton call without running optimisation
newton_guess=zeros(size(waveform));
[x_newton,data]=fmaxnewton(spin_system,@local_quadratic_objective,newton_guess);
result=test_close(result,'fmaxnewton zero iter point',x_newton,newton_guess,0,0,...
                  'fmaxnewton with max_iter zero must return the supplied point');
result=test_true(result,'fmaxnewton zero iter count',...
                 (data.count.iter==0)&&(data.count.gfx==0)&&(data.count.hfx==0),...
                 'fmaxnewton with max_iter zero must not take iterations or derivative calls');

% Check Liouville-space GRAPE derivatives in a tiny vector-space fixture
[liouv_system,drift,controls,rho_init,rho_targ]=local_liouv_fixture('zeeman-liouv');
liouv_wave=[0.7 -0.2; 0.3 0.5];
[~,fid_liouv,grad_liouv,hess_liouv]=grape_liouv(liouv_system,{drift},controls,...
                                                liouv_wave,rho_init,rho_targ,'real');
grad_liouv_ref=local_grape_liouv_grad(liouv_system,drift,controls,...
                                      liouv_wave,rho_init,rho_targ);
result=test_true(result,'grape_liouv finite',...
                 all(isfinite([fid_liouv; grad_liouv(:); hess_liouv(:)])),...
                 'grape_liouv must return finite fidelity, gradient, and Hessian values');
result=test_close(result,'grape_liouv gradient',grad_liouv,grad_liouv_ref,1e-6,1e-6,...
                  'grape_liouv gradient must match centred finite differences');
result=test_close(result,'grape_liouv Hessian symmetry',hess_liouv,hess_liouv.',1e-10,1e-10,...
                  'grape_liouv Hessian must be symmetric for a real fidelity');

% Check TGRAPE duration gradients against finite differences
dt_grid=[0.03; 0.04];
[fid_tgrape,grad_tgrape]=tgrape(liouv_system,drift,controls,liouv_wave,...
                                dt_grid,1,rho_init,rho_targ);
grad_tgrape_ref=local_tgrape_grad(liouv_system,drift,controls,liouv_wave,...
                                  dt_grid,rho_init,rho_targ);
result=test_true(result,'tgrape finite',all(isfinite([fid_tgrape; grad_tgrape(:)])),...
                 'tgrape must return finite fidelity and duration gradient values');
result=test_close(result,'tgrape gradient',grad_tgrape,grad_tgrape_ref,1e-6,1e-6,...
                  'tgrape duration gradients must match centred finite differences');

% Check cooperative phase wrapper on a tiny spherical-tensor-like fixture
coop_system=local_sphten_control_system();
phi_pair=[0.2 -0.3; -0.1 0.4];
[traj_coop,fid_coop,grad_coop]=grape_coop(phi_pair,coop_system); %#ok<ASGLU>
result=test_true(result,'grape_coop finite',all(isfinite([fid_coop(:); grad_coop(:)])),...
                 'grape_coop must return finite cooperative fidelity and gradient channels');
result=test_true(result,'grape_coop gradient shape',isequal(size(grad_coop),[2 2 2]),...
                 'grape_coop must return one phase-gradient row for each cooperative pulse');

end


function local_ensure_pool()

% Start a one-worker process pool if no pool exists
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',1);
end

end


function result=local_check_jacobian(result,label,fun,waveform,J,tol)

% Compare a Jacobian-vector product with centred finite differences
x=waveform(:);
probe=sin((1:numel(x))');
probe=probe/norm(probe,2);
step_size=1e-6;
y_plus=fun(x+step_size*probe);
y_minus=fun(x-step_size*probe);
jv_ref=(y_plus(:)-y_minus(:))/(2*step_size);
jv=J*probe;
result=test_close(result,label,jv,jv_ref,tol,tol,...
                  'the returned Jacobian must predict a finite-difference perturbation');

end


function wave_ref=local_firf_ref(waveform,ker)

% Apply truncated complex convolution one XY channel pair at a time
wave_ref=zeros(size(waveform));
for n=1:(size(waveform,1)/2)
    signal=waveform(2*n-1,:)+1i*waveform(2*n,:);
    filtered=conv(signal(:),ker(:));
    filtered=filtered(1:size(waveform,2));
    wave_ref(2*n-1,:)=real(filtered);
    wave_ref(2*n,:)=imag(filtered);
end

end


function wave_ref=local_spf_ref(waveform,pole)

% Apply the single-pole recurrence explicitly
wave_ref=zeros(size(waveform));
for n=1:(size(waveform,1)/2)
    signal=waveform(2*n-1,:)+1i*waveform(2*n,:);
    filtered=zeros(size(signal));
    filtered(1)=signal(1);
    for k=2:numel(signal)
        filtered(k)=(1-pole)*signal(k)+pole*filtered(k-1);
    end
    wave_ref(2*n-1,:)=real(filtered);
    wave_ref(2*n,:)=imag(filtered);
end

end


function wave_ref=local_szf_ref(waveform,zero)

% Apply the single-zero recurrence explicitly
wave_ref=zeros(size(waveform));
for n=1:(size(waveform,1)/2)
    signal=waveform(2*n-1,:)+1i*waveform(2*n,:);
    filtered=zeros(size(signal));
    filtered(1)=signal(1);
    for k=2:numel(signal)
        filtered(k)=signal(k)/(1-zero)-zero*signal(k-1)/(1-zero);
    end
    wave_ref(2*n-1,:)=real(filtered);
    wave_ref(2*n,:)=imag(filtered);
end

end


function wave_ref=local_amp_tanh_ref(waveform,sat_lvl)

% Apply radial tanh compression while preserving phase
amp=sqrt(waveform(1,:).^2+waveform(2,:).^2);
phi=atan2(waveform(2,:),waveform(1,:));
amp=sat_lvl*tanh(amp/sat_lvl);
wave_ref=[amp.*cos(phi); amp.*sin(phi)];

end


function wave_ref=local_amp_root_ref(waveform,sat_lvl,shape)

% Apply radial root-sigmoid compression while preserving phase
amp=sqrt(waveform(1,:).^2+waveform(2,:).^2);
phi=atan2(waveform(2,:),waveform(1,:));
amp=amp./(1+(amp/sat_lvl).^shape).^(1/shape);
wave_ref=[amp.*cos(phi); amp.*sin(phi)];

end


function wave_ref=local_fapt_ref(fapt,time_grid)

% Add all frequency-amplitude-phase-time events explicitly
wave_ref=zeros(2,numel(time_grid));
for n=1:numel(fapt)
    freq=fapt{n}(1);
    ampl=fapt{n}(2);
    phase=fapt{n}(3);
    start_time=fapt{n}(4);
    end_time=fapt{n}(5);
    time_mask=(time_grid>=start_time)&(time_grid<=end_time);
    wave_ref(1,:)=wave_ref(1,:)+ampl*cos(2*pi*freq*time_grid+phase).*time_mask;
    wave_ref(2,:)=wave_ref(2,:)-ampl*sin(2*pi*freq*time_grid+phase).*time_mask;
end

end


function systems=local_context(~,~,~,~)

% Return two simple context systems for drifts()
H=sparse(diag(1:4));
R=sparse(H);
K=2*sparse(H);
V=3*sparse(H);
systems={{H,R,K},{2*H,R,K,[],V}};

end


function G=local_trap_generator(drifts,controls,dt,c_left,c_right)

% Build the trapezium-product generator explicitly
G_left=drifts{1};
G_right=drifts{2};
for n=1:numel(controls)
    G_left=G_left+c_left(n)*controls{n};
    G_right=G_right+c_right(n)*controls{n};
end
G=(G_left+G_right)/2+1i*dt*(G_left*G_right-G_right*G_left)/12;

end


function D=local_trap_dir(drifts,controls,cc_comm,cc_comm_idx,dt,c_left,c_right,k,side)

% Build a trapezium-product directional derivative explicitly
switch side
    case 'left'
        D=(1/2)*controls{k}+1i*dt*(controls{k}*drifts{2}-drifts{2}*controls{k})/12;
        for n=1:numel(controls)
            if ~cc_comm_idx(n,k)
                D=D+1i*dt*c_right(n)*cc_comm{k,n}/12;
            end
        end
    case 'right'
        D=(1/2)*controls{k}+1i*dt*(drifts{1}*controls{k}-controls{k}*drifts{1})/12;
        for n=1:numel(controls)
            if ~cc_comm_idx(n,k)
                D=D+1i*dt*c_left(n)*cc_comm{n,k}/12;
            end
        end
end

end


function spin_system=local_spin_system(formalism)

% Build a minimal quiet Spinach object
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.sys.scratch=tempdir;
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism=formalism;

end


function spin_system=local_hilb_control_system()

% Build a minimal Hilbert-space one-spin optimal-control system
spin_system=local_spin_system('zeeman-hilb');
S=pauli(2);
control.operators={S.x,S.y};
control.rho_init={S.x};
control.rho_targ={S.z};
control.pwr_levels=1;
control.pulse_dt=[0.02 0.03 0.04];
control.drifts={{sparse(2,2)}};
control.method='lbfgs';
control.max_iter=0;
control.penalties={'none'};
control.p_weights=0;
control.l_bound=-100;
control.u_bound=100;
control.plotting={};
control.parallel='time';
control.amplitudes=[5 6 7];
spin_system=optimcon(spin_system,control);

end


function [spin_system,drift,controls,rho_init,rho_targ]=local_liouv_fixture(formalism)

% Build a tiny vector-space optimal-control fixture
spin_system=local_spin_system(formalism);
spin_system.control.pulse_dt=[0.03 0.04];
spin_system.control.pulse_ntpts=2;
spin_system.control.pulse_nsteps=2;
spin_system.control.integrator='rectangle';
spin_system.control.dead_time=0;
spin_system.control.prefix=[];
spin_system.control.suffix=[];
spin_system.control.keyholes=cell(1,2);
spin_system.control.method='newton';
spin_system.control.plotting={};
spin_system.control.parallel='time';
drift=[0.20 0.10; -0.05 -0.10];
controls={[0 1; 1 0],[0 -1i; 1i 0]};
rho_init=[1; 0.2];
rho_targ=[0.3; 1.0];

end


function spin_system=local_sphten_control_system()

% Build a tiny spherical-tensor-like cooperative-control fixture
spin_system=local_spin_system('sphten-liouv');
spin_system.bas.basis=[0; 1];
control.operators={[0 1; 1 0],[0 -1i; 1i 0]};
control.rho_init={[1; 0.2]};
control.rho_targ={[0.3; 1.0]};
control.pwr_levels=1;
control.pulse_dt=[0.02 0.03];
control.drifts={{sparse(2,2)}};
control.method='lbfgs';
control.max_iter=0;
control.penalties={'none'};
control.p_weights=0;
control.l_bound=-100;
control.u_bound=100;
control.plotting={'correlation_order'};
control.parallel='time';
control.amplitudes=[4 5];
spin_system=optimcon(spin_system,control);

end


function x=local_u2x(u)

% Identity curvilinear-to-Cartesian map
x=u;

end


function J=local_dx_du(~)

% Identity curvilinear-coordinate Jacobian
J=eye(2);

end


function grad=local_phase_grad(phi_profile,spin_system)

% Compute a centred finite-difference phase gradient
step_size=1e-6;
grad=zeros(size(phi_profile));
for n=1:numel(phi_profile)
    phi_plus=phi_profile;
    phi_minus=phi_profile;
    phi_plus(n)=phi_plus(n)+step_size;
    phi_minus(n)=phi_minus(n)-step_size;
    grad(n)=(local_phase_fidelity(phi_plus,spin_system)-...
             local_phase_fidelity(phi_minus,spin_system))/(2*step_size);
end

end


function fidelity=local_phase_fidelity(phi_profile,spin_system)

% Return the physical fidelity channel from grape_phase
[~,fid]=grape_phase(phi_profile,spin_system);
fidelity=fid(1);

end


function [traj_data,fidelity,grad,hess]=local_quadratic_objective(waveform,~)

% Evaluate a simple concave quadratic objective with one dummy penalty
ref_wave=0.1*ones(size(waveform));
residual=waveform-ref_wave;
fidelity=[-sum(residual(:).^2) 0];
traj_data=struct();
if nargout>2
    grad=zeros([size(waveform) 2]);
    grad(:,:,1)=-2*residual;
end
if nargout>3
    hess=zeros(numel(waveform),numel(waveform),2);
    hess(:,:,1)=-2*eye(numel(waveform));
end

end


function grad=local_grape_liouv_grad(spin_system,drift,controls,waveform,rho_init,rho_targ)

% Compute a centred finite-difference Liouville-GRAPE gradient
step_size=1e-6;
grad=zeros(size(waveform));
for n=1:numel(waveform)
    wave_plus=waveform;
    wave_minus=waveform;
    wave_plus(n)=wave_plus(n)+step_size;
    wave_minus(n)=wave_minus(n)-step_size;
    grad(n)=(local_grape_liouv_fid(spin_system,drift,controls,wave_plus,rho_init,rho_targ)-...
             local_grape_liouv_fid(spin_system,drift,controls,wave_minus,rho_init,rho_targ))/(2*step_size);
end

end


function fidelity=local_grape_liouv_fid(spin_system,drift,controls,waveform,rho_init,rho_targ)

% Return the real-overlap fidelity from grape_liouv
[~,fidelity]=grape_liouv(spin_system,{drift},controls,waveform,rho_init,rho_targ,'real');

end


function grad=local_tgrape_grad(spin_system,drift,controls,waveform,dt_grid,rho_init,rho_targ)

% Compute a centred finite-difference TGRAPE duration gradient
step_size=1e-6;
grad=zeros(size(dt_grid));
for n=1:numel(dt_grid)
    dt_plus=dt_grid;
    dt_minus=dt_grid;
    dt_plus(n)=dt_plus(n)+step_size;
    dt_minus(n)=dt_minus(n)-step_size;
    grad(n)=(local_tgrape_fid(spin_system,drift,controls,waveform,dt_plus,rho_init,rho_targ)-...
             local_tgrape_fid(spin_system,drift,controls,waveform,dt_minus,rho_init,rho_targ))/(2*step_size);
end

end


function fidelity=local_tgrape_fid(spin_system,drift,controls,waveform,dt_grid,rho_init,rho_targ)

% Return the TGRAPE fidelity for finite differences
fidelity=tgrape(spin_system,drift,controls,waveform,dt_grid,1,rho_init,rho_targ);

end


