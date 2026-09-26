% Tests frozen input derivatives through waveform maps. Syntax:
%
%                    result=test_freeze_map()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% Noncommuting one-spin transfers exercise temporal filters, phase cycles,
% power scaling, and constrained Hessians in both density formalisms.
%
% talos@spindynamics.org

function result=test_freeze_map()

% Initialise the regression result and a quiet numerical problem
result=new_test_result('optimcon/freeze_map','Frozen waveform-map derivatives',...
                       'Free derivatives must include every physical waveform coordinate.');
spin_system.sys.output='hush'; spin_system.sys.enable={}; spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14; spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5; spin_system.tols.prop_chop=1e-14;
spin_system.comp.isotopes={'1H'}; spin_ops=pauli(2);
formalisms={'zeeman-liouv','zeeman-hilb'};

% Exercise two non-collinear transfers and noncommuting pulse sequences
for fixture=1:2
    if fixture==1
        rho_init=spin_ops.x; rho_targ=spin_ops.y+0.3*spin_ops.z;
        waveform=[3 -1 2;2 4 -3];
    else
        rho_init=spin_ops.y+0.2*spin_ops.z; rho_targ=spin_ops.z-0.4*spin_ops.x;
        waveform=[-4 2 1;1 -3 5];
    end
    for form_idx=1:2

        % Configure equivalent Hilbert and vectorised Liouville problems
        spin_system.bas.formalism=formalisms{form_idx}; control=struct();
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
        control.operators={lx,ly}; control.drifts={{0.7*lz}};
        control.pwr_levels=1.3; control.pulse_dt=[0.04 0.05 0.06];
        control.method='lbfgs'; control.max_iter=0; control.plotting={};
        control.penalties={'none'}; control.p_weights=0;
        control.l_bound=-100; control.u_bound=100;

        % Cover identity, temporal, phase, combined, empty-mask, and trapezium paths
        for variant=1:6
            local_control=control; wave=waveform;
            mask=false(size(wave)); mask(:,end)=true;
            if ismember(variant,[2 4 5 6])
                local_control.distortion={@(w)firf(w,[0.8 0.4]),@(w)spf(w,0.3)};
            end
            if ismember(variant,[3 4 5 6])
                local_control.phase_cycle=[0 0.7 0]; mask(1,1)=true;
            end
            if variant==5, mask(:)=false; end
            if variant==6
                local_control.integrator='trapezium';
                local_control.pulse_dt=[0.04 0.05];
            end
            local_control.freeze=mask;
            local_system=optimcon(spin_system,local_control);
            if variant==5, local_system.control.freeze=[]; end
            [~,fidelity,gradient]=grape_xy(wave,local_system); gradient=gradient(:,:,1);
            label=sprintf('f%d r%d v%d',fixture,form_idx,variant);
            result=test_close(result,[label ' frozen'],gradient(mask),zeros(nnz(mask),1),0,0,...
                              'frozen input derivatives must be exactly zero');

            % Check objective preservation and unchanged unfrozen input derivatives
            free_system=local_system; free_system.control.freeze=[];
            [~,free_fid,free_grad]=grape_xy(wave,free_system); free_grad=free_grad(:,:,1);
            result=test_close(result,[label ' value'],fidelity,free_fid,0,0,...
                              'freezing changes derivatives, not the propagated objective');
            result=test_close(result,[label ' pullback'],gradient(~mask),free_grad(~mask),1e-12,1e-12,...
                              'masking must follow the complete waveform-map pullback');

            % Differentiate the public objective on free input coordinates
            for increment=[1e-3 1e-4]
                numerical=zeros(size(wave));
                for n=find(~mask(:))'
                    plus=wave; minus=wave;
                    plus(n)=plus(n)+increment; minus(n)=minus(n)-increment;
                    [~,fp]=grape_xy(plus,local_system); [~,fm]=grape_xy(minus,local_system);
                    numerical(n)=(fp(1)-fm(1))/(2*increment);
                end
                result=test_close(result,[label ' gradient'],gradient(~mask),numerical(~mask),1e-8,1e-8,...
                                  'active gradients must agree with centred objective differences');
            end
        end

        % Check curvilinear masks after coupled and dimension-changing maps
        for variant=1:9
            local_control=control; wave=waveform;
            u2x=@(u)u; dx_du=@(u)eye(numel(u));
            if variant==2
                u2x=@(u)[u(1)+0.4*u(2);0.3*u(1)+u(2)];
                dx_du=@(u)[1 0.3;0.4 1];
            elseif variant==3
                wave=[waveform;0.5 -0.7 0.2];
                u2x=@(u)[u(1)+0.4*u(3);u(2)+0.3*u(3)];
                dx_du=@(u)[1 0;0 1;0.4 0.3];
            elseif variant==4
                wave=wave(1,:);
                u2x=@(u)[u(1);0.6*u(1)]; dx_du=@(u)[1 0.6];
            elseif variant>=5
                wave=[abs(wave(1,:));0.2*wave(2,:)];
                u2x=@(u)[u(1)*cos(u(2));u(1)*sin(u(2))];
                dx_du=@(u)[cos(u(2)) sin(u(2));-u(1)*sin(u(2)) u(1)*cos(u(2))];
            end
            mask=false(size(wave)); mask(1,1)=true; mask(end,end)=true;
            if ismember(variant,[6 7]), mask(:)=false; end
            if variant==8
                local_control.integrator='trapezium';
                local_control.pulse_dt=[0.04 0.05];
            end
            if variant==9
                local_control.distortion={@(w)firf(w,[0.8 0.4]),@(w)spf(w,0.3)};
                local_control.phase_cycle=[0 0.7 0];
            end
            local_control.penalties={'NS'}; local_control.p_weights=0.2;
            local_control.freeze=mask; local_system=optimcon(spin_system,local_control);
            if variant==6, local_system.control.freeze=[]; end
            label=sprintf('curv f%d r%d v%d',fixture,form_idx,variant);

            % Keep shape failures visible while allowing the remaining fixtures to run
            try
                [~,fidelity,gradient]=grape_curv(wave,u2x,dx_du,local_system);
            catch exception
                result=test_true(result,[label ' call'],false,exception.message);
                continue
            end

            % Compare all objective channels with the unconstrained pullback
            free_system=local_system; free_system.control.freeze=[];
            [~,free_fid,free_grad]=grape_curv(wave,u2x,dx_du,free_system);
            result=test_close(result,[label ' value'],fidelity,free_fid,0,0,...
                              'freezing must leave physical and penalty values unchanged');
            full_mask=repmat(mask,1,1,size(gradient,3));
            result=test_close(result,[label ' frozen'],gradient(full_mask),zeros(nnz(full_mask),1),0,0,...
                              'curvilinear frozen coordinates must vanish in every objective channel');
            result=test_close(result,[label ' pullback'],gradient(~full_mask),free_grad(~full_mask),1e-12,1e-12,...
                              'free curvilinear coordinates must retain all Cartesian contributions');

            % Differentiate the full public objective at two finite-difference increments
            for increment=[1e-4 1e-5]
                numerical=zeros(size(gradient));
                for n=find(~mask(:))'
                    plus=wave; minus=wave;
                    plus(n)=plus(n)+increment; minus(n)=minus(n)-increment;
                    [~,fp]=grape_curv(plus,u2x,dx_du,local_system);
                    [~,fm]=grape_curv(minus,u2x,dx_du,local_system);
                    for k=1:numel(fidelity)
                        numerical(n+(k-1)*numel(wave))=(fp(k)-fm(k))/(2*increment);
                    end
                end
                result=test_close(result,[label ' gradient'],gradient(~full_mask),numerical(~full_mask),1e-8,1e-8,...
                                  'curvilinear derivatives must agree with objective differences');
                fprintf('%s h=%.1e error=%.6e reference=%.6e\n',label,increment,...
                        norm(gradient(~full_mask)-numerical(~full_mask)),norm(numerical(~full_mask)));
            end
        end

        % Preserve phase-only gradients and exact Hessians with power scaling
        local_control=control; local_control.method='newton';
        local_control.amplitudes=[2 3 4]; local_control.freeze=[true false false];
        local_control.phase_cycle=[0 0.7 0];
        local_system=optimcon(spin_system,local_control); phases=[0.2 -0.4 0.7];
        [~,fidelity,gradient,hessian]=grape_phase(phases,local_system);
        free_system=local_system; free_system.control.freeze=[];
        [~,free_fid,free_grad,free_hess]=grape_phase(phases,free_system);
        free_grad(1)=0; free_hess(1,:,1)=0; free_hess(:,1,1)=0;
        result=test_close(result,'phase value',fidelity,free_fid,0,0,'phase freezing preserves values');
        result=test_close(result,'phase gradient',gradient,free_grad,1e-12,1e-12,...
                          'phase freezing preserves free derivatives');
        result=test_close(result,'phase Hessian',hessian,free_hess,1e-12,1e-12,...
                          'phase freezing preserves free curvature');

        % Cover supported exact Hessians with phase mixing and asymmetric freezing
        methods={'newton','goodwin','newton'};
        for method_idx=1:3
            if (form_idx==2)&&(method_idx==3), continue; end
            local_control=control; local_control.method=methods{method_idx};
            local_control.phase_cycle=[0 0.7 0];
            mask=false(size(waveform)); mask(1,1)=true; mask(2,3)=true;
            local_control.freeze=mask;
            if method_idx==3
                local_control.drifts={{0.7*lz-1i*0.3*diag([0 1 1 0])}};
            end
            local_system=optimcon(spin_system,local_control);
            [~,fidelity,gradient,hessian]=grape_xy(waveform,local_system);
            gradient=gradient(:,:,1); hessian=hessian(:,:,1);
            label=sprintf('f%d r%d m%d',fixture,form_idx,method_idx);
            result=test_close(result,[label ' frozen gradient'],gradient(mask),zeros(nnz(mask),1),0,0,...
                              'exact-Hessian calls must preserve frozen gradient zeros');
            result=test_close(result,[label ' frozen rows'],hessian(mask(:),:),zeros(nnz(mask),numel(mask)),0,0,...
                              'frozen input Hessian rows must be zero');
            result=test_close(result,[label ' frozen columns'],hessian(:,mask(:)),zeros(numel(mask),nnz(mask)),0,0,...
                              'frozen input Hessian columns must be zero');
            free_system=local_system; free_system.control.freeze=[];
            [~,free_fid]=grape_xy(waveform,free_system);
            result=test_close(result,[label ' value'],fidelity,free_fid,0,0,...
                              'Hessian masking must not change the objective');

            % Compare active curvature with finite differences of production gradients
            for increment=[1e-3 1e-4]
                numerical=zeros(numel(waveform),numel(waveform));
                for n=find(~mask(:))'
                    plus=waveform; minus=waveform;
                    plus(n)=plus(n)+increment; minus(n)=minus(n)-increment;
                    [~,~,gp]=grape_xy(plus,local_system); [~,~,gm]=grape_xy(minus,local_system);
                    gp=gp(:,:,1); gm=gm(:,:,1); numerical(:,n)=(gp(:)-gm(:))/(2*increment);
                end
                result=test_close(result,[label ' Hessian'],hessian,numerical,1e-8,1e-8,...
                                  'active Hessians must differentiate the constrained gradient');
            end
        end

        % Preserve the direct Liouville engine mask without ensemble maps
        if form_idx==1
            direct_control=control; direct_control.freeze=mask;
            direct_system=optimcon(spin_system,direct_control);
            [~,~,direct_grad]=grape_liouv(direct_system,control.drifts{1},control.operators,...
                                         waveform,control.rho_init{1},control.rho_targ{1},'real');
            direct_system.control.freeze=[];
            [~,~,free_grad]=grape_liouv(direct_system,control.drifts{1},control.operators,...
                                       waveform,control.rho_init{1},control.rho_targ{1},'real');
            free_grad(mask)=0;
            result=test_close(result,'direct engine mask',direct_grad,free_grad,0,0,...
                              'direct engine calls must retain their existing mask semantics');
        else
            direct_control=control; direct_control.freeze=mask;
            direct_system=optimcon(spin_system,direct_control);
            [~,~,direct_grad]=grape_hilb(direct_system,control.drifts{1},control.operators,...
                                        waveform,control.rho_init{1},control.rho_targ{1},'real');
            direct_system.control.freeze=[];
            [~,~,free_grad]=grape_hilb(direct_system,control.drifts{1},control.operators,...
                                      waveform,control.rho_init{1},control.rho_targ{1},'real');
            result=test_close(result,'direct Hilbert mask',direct_grad,free_grad,0,0,...
                              'direct Hilbert calls retain their existing unmasked derivatives');
        end
    end
end

end


