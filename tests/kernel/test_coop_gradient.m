% Tests cooperative phase gradients and the primary fidelity selector.
% Syntax:
%
%                     result=test_coop_gradient()
%
% Outputs:
%
%     result  - regression result with explanatory messages
%
% Two noncommuting pulses are checked in Hilbert and Liouville space,
% with unit and nonunit targets, a complex detection operator, and a
% power ensemble, zero impurity, purely imaginary auxiliary overlaps,
% and cancellation of primary transfer by the impurity penalty, including
% an anonymous forwarding adapter.
% Independent matrix propagation checks the objective;
% centred differences at three increments check its phase gradient.
% All four optimiser methods must reject unusable assembled initial
% guesses while permitting objective-only calls and valid optimisation.
%
% talos@spindynamics.org

function result=test_coop_gradient()

% Initialise the regression result
result=new_test_result('kernel/coop_gradient','Cooperative phase gradients',...
                       'The primary fidelity and squared impurity must share one gradient.');

% Set the small-system numerical environment
spin_system.sys.output='hush';
spin_system.sys.enable={}; spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.comp.isotopes={'1H'};
spin_ops=pauli(2);
formalisms={'zeeman-liouv','zeeman-hilb'};
measures={'real','square'};
methods={'lbfgs','rbfgs','newton','goodwin'};
steps=[1e-3 1e-4 1e-5];

% Exercise a unit target and a nonunit complex detection operator
for fixture=1:4
    amplitudes=[3 4 2]; pulse_dt=[0.04 0.05 0.06]; drift_scale=0.7;
    if fixture==1
        rho_init=spin_ops.x;
        rho_targ=spin_ops.x+0.3*spin_ops.z;
        rho_targ=rho_targ/norm(rho_targ,'fro');
        phase_pair=[0.2 -0.3 0.5;-0.1 0.4 0.7];
        powers=1;
    elseif fixture==2
        rho_init=spin_ops.y+0.2*spin_ops.z;
        rho_targ=1.7*(spin_ops.z-0.4*spin_ops.x+1i*spin_ops.y);
        phase_pair=[-0.6 0.9 -0.2;0.8 -0.5 1.1];
        powers=[0.8 1.1];
    elseif fixture==3
        rho_init=[0 1;0 0]; rho_targ=[1 3;1+1i 0];
        phase_pair=[0 pi;0 0]; powers=1;
        amplitudes=[1 1]; pulse_dt=[pi/2 pi/2]; drift_scale=0;
    else
        rho_init=eye(2); rho_targ=eye(2);
        phase_pair=[0.2 -0.3;-0.1 0.4]; powers=1;
        amplitudes=[0 0]; pulse_dt=[0.1 0.2]; drift_scale=0;
    end

    % Independently propagate both experiments in Hilbert space
    overlaps=zeros(numel(powers),2); dirt_cost=zeros(numel(powers),1);
    for power_idx=1:numel(powers)
        dirt_sum=zeros(2);
        for pulse_idx=1:2
            rho=rho_init;
            for slice=1:numel(pulse_dt)
                H=drift_scale*spin_ops.z+powers(power_idx)*amplitudes(slice)*...
                  (cos(phase_pair(pulse_idx,slice))*spin_ops.x+...
                   sin(phase_pair(pulse_idx,slice))*spin_ops.y);
                P=expm(-1i*H*pulse_dt(slice)); rho=P*rho*P';
            end
            overlap=sum(conj(rho_targ).*rho,'all');
            overlaps(power_idx,pulse_idx)=overlap;
            dirt_sum=dirt_sum+rho-rho_targ*overlap/norm(rho_targ,'fro')^2;
        end
        dirt_cost(power_idx)=norm(dirt_sum,'fro')^2;
    end

    % Configure equivalent density-operator representations
    for form_idx=1:2
        spin_system.bas.formalism=formalisms{form_idx};
        control=struct();
        control.isotopes={'1H'}; control.channels=[1;1];
        if form_idx==1
            lx=kron(eye(2),spin_ops.x)-kron(transpose(spin_ops.x),eye(2));
            ly=kron(eye(2),spin_ops.y)-kron(transpose(spin_ops.y),eye(2));
            lz=kron(eye(2),spin_ops.z)-kron(transpose(spin_ops.z),eye(2));
            control.rho_init={rho_init(:)}; control.rho_targ={rho_targ(:)};
        else
            lx=spin_ops.x; ly=spin_ops.y; lz=spin_ops.z;
            control.rho_init={rho_init}; control.rho_targ={rho_targ};
        end
        control.operators={lx,ly}; control.drifts={{drift_scale*lz}};
        control.pwr_levels=powers; control.pulse_dt=pulse_dt;
        control.method='lbfgs'; control.max_iter=0;
        control.penalties={'none'}; control.p_weights=0;
        control.l_bound=-100; control.u_bound=100; control.plotting={};
        control.amplitudes=amplitudes;

        % Preserve each primary fidelity while penalising squared impurity
        for measure_idx=1:2
            control.fidelity=measures{measure_idx};
            local_system=optimcon(spin_system,control);
            [traj_data,fidelity,gradient]=grape_coop(phase_pair,local_system);
            gradient=gradient(:,:,1);
            switch measures{measure_idx}
                case 'real', primary=real(overlaps);
                case 'square', primary=abs(overlaps).^2;
            end
            expected=mean(primary,'all')-mean(dirt_cost);
            label=sprintf('fixture %d %s %s',fixture,formalisms{form_idx},measures{measure_idx});
            result=test_close(result,[label ' objective'],fidelity(1),expected,1e-11,0,...
                              'Independent propagation must retain the requested primary measure.');

            % Keep both trajectory branches and report the pre-impurity primary
            [result,has_primary]=test_true(result,[label ' trajectory metadata'],...
                                           numel(traj_data)==2&&iscell(traj_data{1})&&...
                                           ~isempty(traj_data{1})&&isstruct(traj_data{1}{1})&&...
                                           isfield(traj_data{1}{1},'primary_fid'),...
                                           'The first trajectory must carry the primary transfer.');
            if has_primary
                result=test_close(result,[label ' primary transfer'],...
                                  traj_data{1}{1}.primary_fid,mean(primary,'all'),1e-11,0,...
                                  'The optimiser must receive transfer before impurity subtraction.');
            end

            % Compare every phase derivative at three finite-difference increments
            for step_size=steps
                fd_grad=zeros(size(phase_pair));
                for n=1:numel(phase_pair)
                    plus=phase_pair; minus=phase_pair;
                    plus(n)=plus(n)+step_size; minus(n)=minus(n)-step_size;
                    [~,fp]=grape_coop(plus,local_system);
                    [~,fm]=grape_coop(minus,local_system);
                    fd_grad(n)=(fp(1)-fm(1))/(2*step_size);
                end
                fprintf('COOP %s h=%.1e error=%.12e scale=%.12e\n',...
                        label,step_size,norm(gradient-fd_grad,'fro'),norm(fd_grad,'fro'));

                % Scale the second-order error allowance with the reference gradient
                result=test_close(result,sprintf('%s h=%.1e',label,step_size),...
                                  gradient,fd_grad,2e-9,2*step_size^2,...
                                  'The derivative must match centred differences to second order.');
            end
        end

        % Check exact zeros without confusing a costate with a primary target
        if fixture==3
            rho_a=[0 1;0 0]; rho_b=[0 0;1 0];
            dirt_sum=rho_a+rho_b-rho_targ*hdot(rho_targ,rho_a+rho_b)/hdot(rho_targ,rho_targ);
            overlap=hdot(dirt_sum,rho_a);
            result=test_true(result,[formalisms{form_idx} ' imaginary auxiliary'],...
                             real(overlap)==0&&imag(overlap)~=0,...
                             'A valid auxiliary overlap is exactly nonzero and purely imaginary.');
            if form_idx==1
                engine=@grape_liouv; source=rho_a(:); target=dirt_sum(:);
                identity=reshape(eye(2),[],1);
            else
                engine=@grape_hilb; source=rho_a; target=dirt_sum; identity=eye(2);
            end
            waveform=zeros(2,2);
            [~,aux_fid,aux_grad]=engine(local_system,control.drifts{1},...
                                      control.operators,waveform,source,target,'real');
            result=test_close(result,[formalisms{form_idx} ' zero real auxiliary'],...
                              aux_fid,0,0,0,'Zero real overlap must return its nonzero derivative.');
            result=test_true(result,[formalisms{form_idx} ' nonzero auxiliary gradient'],...
                             norm(aux_grad,'fro')>0,'The real-linear auxiliary derivative is nonzero.');

            % Difference independent matrix exponentials, not engine fidelities
            for step_size=steps
                fd_grad=zeros(size(waveform));
                for n=1:numel(waveform)
                    values=zeros(1,2);
                    for side=1:2
                        perturbed=waveform; perturbed(n)=(3-2*side)*step_size;
                        rho=rho_a;
                        for slice=1:numel(pulse_dt)
                            H=perturbed(1,slice)*spin_ops.x+perturbed(2,slice)*spin_ops.y;
                            P=expm(-1i*H*pulse_dt(slice)); rho=P*rho*P';
                        end
                        values(side)=real(sum(conj(dirt_sum).*rho,'all'));
                    end
                    fd_grad(n)=(values(1)-values(2))/(2*step_size);
                end
                fprintf('AUX %s h=%.1e error=%.12e scale=%.12e\n',...
                        formalisms{form_idx},step_size,norm(aux_grad-fd_grad,'fro'),norm(fd_grad,'fro'));
                result=test_close(result,sprintf('%s auxiliary h=%.1e',formalisms{form_idx},step_size),...
                                  aux_grad,fd_grad,2*step_size^2+2e-9,0,...
                                  'The zero-value auxiliary derivative must match independent propagation.');
            end

            % Vanishing impurity has a vanishing derivative in both formalisms
            [~,zero_fid,zero_grad]=engine(local_system,control.drifts{1},...
                                        control.operators,waveform,source,zeros(size(target)),'real');
            result=test_close(result,[formalisms{form_idx} ' zero impurity'],...
                              [zero_fid;zero_grad(:)],zeros(5,1),0,0,...
                              'A zero impurity costate must give exact zero value and gradient.');
            [~,const_fid,const_grad]=engine(local_system,control.drifts{1},...
                                          control.operators,waveform,identity,identity,'real');
            result=test_close(result,[formalisms{form_idx} ' constant overlap'],...
                              [const_fid;const_grad(:)],[2;zeros(4,1)],0,0,...
                              'A nonzero constant overlap must have an exact zero gradient.');

            % Check the assembled-objective safeguard for every optimiser method
            local_system.control.fidelity='real'; local_system.control.pulse_dt=[0.1 0.2];
            for method_idx=1:numel(methods)
                local_system.control.method=methods{method_idx};
                local_system.control.max_iter=1; local_system.control.freeze=[];
                label=[formalisms{form_idx} ' ' methods{method_idx}];
                for guard_case=1:2
                    if guard_case==1
                        local_system.control.rho_init={source};
                        local_system.control.rho_targ={target};
                    else
                        local_system.control.rho_init={identity};
                        local_system.control.rho_targ={identity};
                    end
                    caught=''; lastwarn('');
                    try
                        fmaxnewton(local_system,@grape_xy,waveform);
                    catch err
                        caught=err.message;
                    end
                    warn_text=lastwarn;
                    result=test_true(result,sprintf('%s optimiser guard %d',label,guard_case),...
                                     strcmp(caught,'fidelity or gradient too small at iter 1, find a better guess.')&&...
                                     isempty(warn_text),'An unusable initial guess must fail before a singular solve.');
                end

                % Preserve objective-only evaluation of constant objectives
                local_system.control.max_iter=0;
                [point,data]=fmaxnewton(local_system,@grape_xy,waveform);
                result=test_true(result,[label ' objective only'],isequal(point,waveform)&&...
                                 data.count.fx==1&&data.count.gfx==0&&data.count.hfx==0,...
                                 'Zero iterations must not request derivatives or reject a constant objective.');

                % Optimise a nonstationary physical transfer with each supported method
                source_ref=spin_ops.x; target_ref=spin_ops.x+0.3*spin_ops.z;
                if form_idx==1
                    local_system.control.rho_init={source_ref(:)};
                    local_system.control.rho_targ={target_ref(:)};
                else
                    local_system.control.rho_init={source_ref};
                    local_system.control.rho_targ={target_ref};
                end
                local_system.control.max_iter=3; guess=[0.2 -0.3;0.4 0.5];
                [~,before]=grape_xy(guess,local_system);
                [point,data]=fmaxnewton(local_system,@grape_xy,guess);
                [~,after]=grape_xy(point,local_system);
                fprintf('OPTIMISER %s before=%.15g after=%.15g iter=%d hfx=%d\n',...
                        label,before(1),after(1),data.count.iter,data.count.hfx);
                result=test_true(result,[label ' nonzero objective'],all(isfinite(point),'all')&&...
                                 after(1)>before(1),'A valid initial guess must still improve its transfer fidelity.');
                result=test_true(result,[label ' derivative counts'],data.count.iter>1&&...
                                 data.count.hfx==data.count.iter*ismember(methods{method_idx},{'newton','goodwin'}),...
                                 'Only Hessian methods request Hessians, once per iteration.');

                % Refuse a nonzero gradient when every input coordinate is frozen
                local_system.control.freeze=true(size(guess)); caught=''; lastwarn('');
                try
                    fmaxnewton(local_system,@grape_xy,guess);
                catch err
                    caught=err.message;
                end
                warn_text=lastwarn;
                result=test_true(result,[label ' frozen gradient'],...
                                 strcmp(caught,'fidelity or gradient too small at iter 1, find a better guess.')&&...
                                 isempty(warn_text),'Only unfrozen coordinates may contribute to the initial gradient.');
            end
        end
    end
end

% Check a nonstationary cooperative objective at a zero-score cancellation
spin_system.bas.formalism='zeeman-hilb';
coop_control=struct();
coop_control.isotopes={'1H'}; coop_control.channels=[1;1];
coop_control.operators={spin_ops.x,spin_ops.y};
coop_control.drifts={{0.1*spin_ops.z}};
coop_control.rho_init={spin_ops.x};
coop_control.rho_targ={(spin_ops.x+0.3*spin_ops.z)/norm(spin_ops.x+0.3*spin_ops.z,'fro')};
coop_control.pwr_levels=1; coop_control.pulse_dt=[0.2 0.2];
coop_control.method='lbfgs'; coop_control.max_iter=1;
coop_control.penalties={'none'}; coop_control.p_weights=0;
coop_control.l_bound=-100; coop_control.u_bound=100;
coop_control.plotting={}; coop_control.amplitudes=[10 10];
coop_system=optimcon(spin_system,coop_control);
coop_guess=0.364110104613*ones(2,2);
[~,coop_before,coop_gradient]=grape_coop(coop_guess,coop_system);
[~,primary_a]=grape_phase(coop_guess(1,:),coop_system);
[~,primary_b]=grape_phase(coop_guess(2,:),coop_system);
primary=(primary_a(1)+primary_b(1))/2;
result=test_true(result,'cooperative cancellation fixture',...
                 abs(coop_before(1))<1e-6&&primary>1e-3&&...
                 norm(coop_gradient(:,:,1),'fro')>1e-3,...
                 'The nearly zero composite score hides substantial primary transfer and gradient.');
coop_caught='';
try
    [coop_point,coop_data]=fmaxnewton(coop_system,@grape_coop,coop_guess);
catch err
    coop_caught=err.message;
end
result=test_true(result,'cooperative cancellation admitted',isempty(coop_caught),...
                 'The optimiser must not confuse impurity cancellation with zero primary transfer.');
if isempty(coop_caught)
    [~,coop_after]=grape_coop(coop_point,coop_system);
    result=test_true(result,'cooperative cancellation improves',...
                     coop_after(1)>coop_before(1)&&coop_data.count.iter==1,...
                     'A nonstationary zero-score waveform must improve in a real optimiser step.');
end

% Check identical cooperative admission through an anonymous adapter
coop_adapter=@(wave,system)grape_coop(wave,system);
coop_caught='';
try
    [adapt_point,adapt_data]=fmaxnewton(coop_system,coop_adapter,coop_guess);
catch err
    coop_caught=err.message;
end
result=test_true(result,'forwarded cooperative cancellation admitted',isempty(coop_caught),...
                 'Callback forwarding must not turn a composite score into primary transfer.');
if isempty(coop_caught)
    [~,adapt_after]=grape_coop(adapt_point,coop_system);
    result=test_true(result,'forwarded cooperative cancellation improves',...
                     adapt_after(1)>coop_before(1)&&adapt_data.count.iter==1&&...
                     adapt_data.count.fx==coop_data.count.fx,...
                     'The forwarded objective must improve in a real optimiser step.');
end

% Reject a nonstationary impurity direction with no primary transfer
zero_system=coop_system; zero_system.control.rho_targ={eye(2)};
zero_guess=coop_guess; zero_guess(2,1)=zero_guess(2,1)+pi/4;
[zero_traj,~,zero_gradient]=grape_coop(zero_guess,zero_system);
result=test_true(result,'zero primary cooperative fixture',...
                 abs(zero_traj{1}{1}.primary_fid)<1e-12&&...
                 norm(zero_gradient(:,:,1),'fro')>1e-3,...
                 'A nonzero impurity gradient must not conceal zero primary transfer.');
zero_caught='';
try
    fmaxnewton(zero_system,@grape_coop,zero_guess);
catch err
    zero_caught=err.message;
end
result=test_true(result,'zero primary cooperative rejected',...
                 strcmp(zero_caught,'fidelity or gradient too small at iter 1, find a better guess.'),...
                 'The optimiser must retain the true primary-fidelity safeguard.');

% Keep the assembled-gradient guard when every cooperative phase is frozen
coop_system.control.freeze=true(size(coop_guess)); coop_caught='';
try
    fmaxnewton(coop_system,@grape_coop,coop_guess);
catch err
    coop_caught=err.message;
end
result=test_true(result,'frozen cooperative gradient rejected',...
                 strcmp(coop_caught,'fidelity or gradient too small at iter 1, find a better guess.'),...
                 'The value exemption must not weaken the assembled-gradient guard.');

end


