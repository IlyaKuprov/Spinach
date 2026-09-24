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
% power ensemble. Independent matrix propagation checks the objective;
% centred differences at three increments check its phase gradient.
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
steps=[1e-3 1e-4 1e-5];

% Exercise a unit target and a nonunit complex detection operator
for fixture=1:2
    if fixture==1
        rho_init=spin_ops.x;
        rho_targ=spin_ops.x+0.3*spin_ops.z;
        rho_targ=rho_targ/norm(rho_targ,'fro');
        phase_pair=[0.2 -0.3 0.5;-0.1 0.4 0.7];
        powers=1;
    else
        rho_init=spin_ops.y+0.2*spin_ops.z;
        rho_targ=1.7*(spin_ops.z-0.4*spin_ops.x+1i*spin_ops.y);
        phase_pair=[-0.6 0.9 -0.2;0.8 -0.5 1.1];
        powers=[0.8 1.1];
    end

    % Independently propagate both experiments in Hilbert space
    amplitudes=[3 4 2]; pulse_dt=[0.04 0.05 0.06];
    overlaps=zeros(numel(powers),2); dirt_cost=zeros(numel(powers),1);
    for power_idx=1:numel(powers)
        dirt_sum=zeros(2);
        for pulse_idx=1:2
            rho=rho_init;
            for slice=1:3
                H=0.7*spin_ops.z+powers(power_idx)*amplitudes(slice)*...
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
        control.operators={lx,ly}; control.drifts={{0.7*lz}};
        control.pwr_levels=powers; control.pulse_dt=pulse_dt;
        control.method='lbfgs'; control.max_iter=0;
        control.penalties={'none'}; control.p_weights=0;
        control.l_bound=-100; control.u_bound=100; control.plotting={};
        control.amplitudes=amplitudes;

        % Preserve each primary fidelity while penalising squared impurity
        for measure_idx=1:2
            control.fidelity=measures{measure_idx};
            local_system=optimcon(spin_system,control);
            [~,fidelity,gradient]=grape_coop(phase_pair,local_system);
            gradient=gradient(:,:,1);
            switch measures{measure_idx}
                case 'real', primary=real(overlaps);
                case 'square', primary=abs(overlaps).^2;
            end
            expected=mean(primary,'all')-mean(dirt_cost);
            label=sprintf('fixture %d %s %s',fixture,formalisms{form_idx},measures{measure_idx});
            result=test_close(result,[label ' objective'],fidelity(1),expected,1e-11,0,...
                              'Independent propagation must retain the requested primary measure.');

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
                result=test_close(result,sprintf('%s h=%.1e',label,step_size),...
                                  gradient,fd_grad,2*step_size^2+2e-9,0,...
                                  'The derivative must match centred differences to second order.');
            end
        end
    end
end

end


