% Tests the supported boundaries of keyhole optimal control. Syntax:
%
%                    result=test_keyhole_guard()
%
% Outputs:
%
%     result  - regression test results and explanatory messages
%
% Spin-state projectors act between noncommuting spin-half pulse slices.
% Density-vector and wavefunction bases reject Newton/Goodwin keyhole
% methods at setup and direct-engine entry, retaining first derivatives.
% Empty schedules and Hilbert-space keyhole Hessians remain supported.
% The spherical basis uses normalised identity, T11, T10, and T1-1.
%
% talos@spindynamics.org

function result=test_keyhole_guard()

% Declare the physically motivated regression
result=new_test_result('kernel/keyhole_guard','State-vector keyhole Hessian guard',...
                       'Unsupported Hessians must be refused without changing supported methods.');

% Define a spin-half system with complex noncommuting controls
spin_system.sys.output='hush';
spin_system.sys.enable={}; spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.comp.isotopes={'1H'};
spin_ops=pauli(2);
rho_init=spin_ops.y+0.2*spin_ops.z;
rho_targ=spin_ops.z-0.4*spin_ops.x;
waveform=[-4 2 1;1 -3 5];
methods={'lbfgs','rbfgs','newton','goodwin'};
formalisms={'zeeman-liouv','sphten-liouv','zeeman-hilb','zeeman-wavef'};

% Change coordinates explicitly into normalised spherical tensors
Q=[1/sqrt(2) 0 1/sqrt(2) 0;...
   0 0 0 1;0 -1 0 0;1/sqrt(2) 0 -1/sqrt(2) 0];

% Exercise the density-operator and wavefunction formalisms
for form_idx=1:numel(formalisms)
    spin_system.bas.formalism=formalisms{form_idx};
    control=struct();
    control.isotopes={'1H'}; control.channels=[1;1];
    if form_idx==3
        controls={spin_ops.x,spin_ops.y}; drift=0.7*spin_ops.z;
        source=rho_init; target=rho_targ;
        keyhole=@(rho)diag(diag(rho));
    elseif form_idx==4
        controls={spin_ops.x,spin_ops.y}; drift=0.7*spin_ops.z;
        source=[2;1]/sqrt(5); target=[1;-1i]/sqrt(2); P=diag([1 0]);
        keyhole=@(rho)P*rho; control.fidelity='square';
    else
        controls={hilb2liouv(spin_ops.x,'comm'),hilb2liouv(spin_ops.y,'comm')};
        drift=0.7*hilb2liouv(spin_ops.z,'comm');
        source=rho_init(:); target=rho_targ(:); P=diag([1 0 0 1]);
        if form_idx==2
            controls={Q'*controls{1}*Q,Q'*controls{2}*Q};
            drift=Q'*drift*Q; source=Q'*source; target=Q'*target; P=Q'*P*Q;
        end
        keyhole=@(rho)P*rho;
    end

    % Supply the complete minimal control problem
    control.operators=controls; control.drifts={{drift}};
    control.rho_init={source}; control.rho_targ={target};
    control.pwr_levels=1; control.pulse_dt=[0.04 0.05 0.06];
    control.max_iter=0; control.penalties={'none'}; control.p_weights=0;
    control.l_bound=-100; control.u_bound=100; control.plotting={};

    % Check first-order projection derivatives and every exact-Hessian boundary
    for method_idx=1:numel(methods)
        control.method=methods{method_idx};
        control.keyholes={keyhole,keyhole,[]};
        label=[formalisms{form_idx} ' ' methods{method_idx}];
        if (form_idx~=3)&&(method_idx>2)

            % Require the specific unsupported-combination error from setup
            error_text='';
            try
                optimcon(spin_system,control);
            catch exception
                error_text=exception.message;
            end
            result=test_true(result,[label ' setup guard'],...
                             contains(error_text,'keyholes with Newton/Goodwin Hessians')&&...
                             contains(error_text,'not implemented'),...
                             'setup must explicitly refuse the unsupported combination');

            % Bypass setup method validation without changing the physical problem
            control.method='lbfgs';
            local_system=optimcon(spin_system,control);
            local_system.control.method=methods{method_idx};
            error_text='';
            try
                [~,~,~,~]=grape_liouv(local_system,{drift},controls,...
                                     waveform,source,target,local_system.control.fidelity);
            catch exception
                error_text=exception.message;
            end
            result=test_true(result,[label ' direct guard'],...
                             contains(error_text,'keyholes with Newton/Goodwin Hessians')&&...
                             contains(error_text,'not implemented'),...
                             'direct engine calls must not bypass the unsupported boundary');

            % Keep empty schedules available for both exact-Hessian methods
            control.method=methods{method_idx}; control.keyholes=cell(1,3);
        end

        % Preserve the requested method and compute supported derivatives
        local_system=optimcon(spin_system,control);
        result=test_true(result,[label ' method retained'],...
                         strcmp(local_system.control.method,methods{method_idx}),...
                         'setup must never substitute another optimisation algorithm');
        if (form_idx~=3)&&(method_idx<=2)

            % Refuse requested Hessians even under first-order methods
            error_text='';
            try
                [~,~,~,~]=grape_xy(waveform,local_system);
            catch exception
                error_text=exception.message;
            end
            result=test_true(result,[label ' requested Hessian guard'],...
                             ~isempty(error_text),...
                             'first-order methods must not expose incorrect keyhole Hessians');

            % Enforce the same boundary on direct engine calls
            error_text='';
            try
                [~,~,~,~]=grape_liouv(local_system,{drift},controls,...
                                     waveform,source,target,local_system.control.fidelity);
            catch exception
                error_text=exception.message;
            end
            result=test_true(result,[label ' direct Hessian guard'],...
                             contains(error_text,'keyhole Hessians')&&...
                             contains(error_text,'not implemented'),...
                             'direct requests must not expose incorrect keyhole Hessians');
        end
        if method_idx>2
            [~,~,gradient,hessian]=grape_xy(waveform,local_system);
        else
            [~,~,gradient]=grape_xy(waveform,local_system);
        end
        gradient=gradient(:,:,1);

        % Independently differentiate the composed physical objective
        step_size=1e-4; grad_ref=zeros(size(waveform));
        hess_ref=zeros(numel(waveform),numel(waveform));
        for n=1:numel(waveform)
            plus=waveform; minus=waveform;
            plus(n)=plus(n)+step_size; minus(n)=minus(n)-step_size;
            [~,fp,gp]=grape_xy(plus,local_system);
            [~,fm,gm]=grape_xy(minus,local_system);
            grad_ref(n)=(fp(1)-fm(1))/(2*step_size);
            gp=gp(:,:,1); gm=gm(:,:,1);
            hess_ref(:,n)=(gp(:)-gm(:))/(2*step_size);
        end
        result=test_close(result,[label ' gradient'],gradient,grad_ref,1e-9,0,...
                          'first-order projected and supported Hessian paths retain their gradients');
        if method_idx>2
            result=test_close(result,[label ' Hessian'],hessian(:,:,1),hess_ref,1e-9,0,...
                              'empty vector schedules and projected Hilbert Hessians remain valid');
        end
    end

    % Refuse dissipative Newton keyholes without relying on unitary dynamics
    if form_idx<3
        control.method='newton'; control.keyholes={keyhole,keyhole,[]};
        control.drifts={{drift-1i*0.3*(eye(4)-P)}};
        error_text='';
        try
            optimcon(spin_system,control);
        catch exception
            error_text=exception.message;
        end
        result=test_true(result,[formalisms{form_idx} ' damped guard'],...
                         contains(error_text,'keyholes with Newton/Goodwin Hessians')&&...
                         contains(error_text,'not implemented'),...
                         'dissipative Newton keyholes are also explicitly unsupported');
    end
end

end


