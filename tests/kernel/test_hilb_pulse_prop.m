% Tests Hilbert shaped-pulse propagators and density-state reuse. Syntax:
%
%                    result=test_hilb_pulse_prop()
%
% Outputs:
%
%     result - regression checks for all six method/quadrature choices
%
% Noncommuting spin-half pulses, constant pulses, and zero durations are
% compared with ordered dense exponentials. Both Hilbert step algorithms
% are exercised; one-sided wavefunction and Liouville calls are controls.
%
% talos@spindynamics.org

function result=test_hilb_pulse_prop()

% Describe the regression target
result=new_test_result('kernel/hilb_pulse_prop',...
                       'Hilbert shaped-pulse propagators',...
                       'Returned propagators must reproduce two-sided density evolution.');

% Build a quiet spin-half system
sys.magnet=0; sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb'; bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Define complex noncommuting generators and two positive density matrices
ops=pauli(2); controls={ops.x,ops.y};
drift=2*pi*(35*ops.z+17*eye(2));
rho=eye(2)/2+0.3*ops.z+0.2*ops.y;
rho_other=eye(2)/2-0.4*ops.x+0.1*ops.y;
methods={'expv-pwc','expv-pwl','expm-pwc','expm-pwl','evol-pwc','evol-pwl'};

% Exercise both small-matrix and commutator-series Hilbert state paths
for cutoff=[100 0]
    spin_system.tols.small_matrix=cutoff;
    for fixture=1:3
        for n=1:numel(methods)

            % Select noncommuting, constant, or zero-duration pulses
            amp_x=2*pi*[430 150 -100 240];
            amp_y=2*pi*[0 220 310 -70];
            slice_durs=[0.0001 0.0002 0.00015];
            if fixture==2
                amp_x(:)=amp_x(1); amp_y(:)=amp_y(1);
            elseif fixture==3
                slice_durs(:)=0;
            end
            if strcmp(methods{n}(6:8),'pwc')
                amp_x=amp_x(1:3); amp_y=amp_y(1:3);
            end

            % Build the ordered reference and all trajectory points
            ref_prop=eye(2); ref_traj=cell(1,4); ref_traj{1}=rho;
            for k=1:numel(slice_durs)
                H=drift+amp_x(k)*ops.x+amp_y(k)*ops.y;
                if strcmp(methods{n}(6:8),'pwl')
                    right=drift+amp_x(k+1)*ops.x+amp_y(k+1)*ops.y;
                    H=(H+right)/2+(1i*slice_durs(k)/6)*(H*right-right*H);
                end
                ref_prop=expm(-1i*H*slice_durs(k))*ref_prop;
                ref_traj{k+1}=ref_prop*rho*ref_prop';
            end

            % Request the state, trajectory, and reusable propagator
            [observed,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...
                                            {amp_x,amp_y},slice_durs,rho,methods{n});
            label=sprintf('%s cutoff=%g fixture=%d',methods{n},cutoff,fixture);
            result=test_close(result,[label ' operator'],P,ref_prop,1e-9,1e-9,...
                              'the returned operator is the ordered one-sided product');
            result=test_close(result,[label ' state'],observed,ref_traj{end},1e-9,1e-9,...
                              'density-state evolution remains two-sided');
            result=test_close(result,[label ' reuse'],P*rho*P',observed,1e-9,1e-9,...
                              'the returned propagator reproduces the returned density');
            result=test_close(result,[label ' unitary'],P'*P,eye(2),1e-9,1e-9,...
                              'Hermitian spin Hamiltonians generate unitary propagators');
            result=test_close(result,[label ' trace'],trace(observed),1,1e-9,1e-9,...
                              'unitary evolution preserves density trace');
            result=test_close(result,[label ' Hermitian'],observed,observed',1e-9,1e-9,...
                              'unitary evolution preserves density Hermiticity');
            for k=1:numel(traj)
                result=test_close(result,[label ' trajectory ' num2str(k)],...
                                  traj{k},ref_traj{k},1e-9,1e-9,...
                                  'every trajectory point retains the two-sided state action');
            end

            % Verify output-count invariance and reuse on a different density
            one_out=shaped_pulse_xy(spin_system,drift,controls,...
                                   {amp_x,amp_y},slice_durs,rho,methods{n});
            [two_out,two_traj]=shaped_pulse_xy(spin_system,drift,controls,...
                                             {amp_x,amp_y},slice_durs,rho,methods{n});
            [other,~,Q]=shaped_pulse_xy(spin_system,drift,controls,...
                                      {amp_x,amp_y},slice_durs,rho_other,methods{n});
            result=test_close(result,[label ' one output'],one_out,observed,0,0,...
                              'requesting a propagator does not change the state computation');
            result=test_close(result,[label ' two outputs'],two_out,two_traj{end},0,0,...
                              'two-output calls retain the final trajectory state');
            result=test_close(result,[label ' output count'],two_out,observed,0,0,...
                              'three-output calls retain the same density-state algorithm');
            result=test_close(result,[label ' other density'],P*rho_other*P',other,1e-9,1e-9,...
                              'the propagator is reusable on a different initial density');
            result=test_close(result,[label ' independent'],P,Q,0,0,...
                              'the propagator does not depend on the initial density');
        end
    end
end

% Preserve the one-sided formalisms across all methods and quadratures
for formalism={'zeeman-wavef','zeeman-liouv'}
    bas.formalism=formalism{1};
    spin_system=test_spin_system(sys,inter,bas);
    controls={operator(spin_system,'Lx',1),operator(spin_system,'Ly',1)};
    drift=2*pi*35*operator(spin_system,'Lz',1);
    initial=ones(size(drift,1),1)/sqrt(size(drift,1));
    for n=1:numel(methods)

        % Use a nonzero constant pulse as the unchanged-formalism control
        slice_durs=[0.0001 0.0002];
        count=2+strcmp(methods{n}(6:8),'pwl');
        amplitudes={2*pi*430*ones(1,count),2*pi*170*ones(1,count)};
        H=drift+2*pi*430*controls{1}+2*pi*170*controls{2};
        ref_prop=expm(full(-1i*H*sum(slice_durs)));
        [observed,~,P]=shaped_pulse_xy(spin_system,drift,controls,...
                                     amplitudes,slice_durs,initial,methods{n});
        label=[formalism{1} ' ' methods{n}];
        result=test_close(result,[label ' operator'],P,ref_prop,1e-8,1e-8,...
                          'one-sided formalisms retain their effective propagator');
        result=test_close(result,[label ' reuse'],P*initial,observed,1e-8,1e-8,...
                          'one-sided formalisms retain state-vector reuse');
    end
end

end


