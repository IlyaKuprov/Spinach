% Comparison of explicitly sliced and Fokker-Planck MAS evolution
% for single-crystal 13C CSA and 27Al quadrupolar systems. The 13C
% calculation includes RF; the 27Al central transition includes the
% third-order rotating-frame correction.
%
% Syntax: mas_fplanck_slices()
%
% Prints convergence metrics and fails if the two routes disagree.
%
% Calculation time: minutes
%
% talos@spindynamics.org

function mas_fplanck_slices()

% 13C system specification
sys.magnet=9.4;
sys.isotopes={'13C'};
sys.parallel={'processes',1};
sys.parprops={};
inter.zeeman.eigs={[-120 -25 145]};
inter.zeeman.euler={[0.4 0.7 0.2]};

% Basis set
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spin system
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% MAS and RF parameters
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='single_crystal';
parameters.spins={'13C'};
parameters.offset=0;
parameters.rframes={{'13C',1}};
parameters.rho0=state(spin_system,'L+','13C');
parameters.coil=state(spin_system,'Lz','13C');
parameters.ref_norm=norm(parameters.rho0)*norm(parameters.coil);
parameters.rf_op=operator(spin_system,'Lx','13C');
parameters.rf_amp=2*pi*3500;
parameters.duration=1/parameters.rate;
parameters.serial=true;
parameters.verbose=0;

% Compare independent rotor-rank and slice-count refinements
[fp_csa,sl_csa]=compare_routes(spin_system,parameters,...
                               [4 6 8 10],[33 65 129]);
check_limit('13C CSA with RF',fp_csa,sl_csa,1e-3);

% 27Al quadrupolar system specification
clear sys inter bas
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'27Al'};
sys.parallel={'processes',1};
sys.parprops={};
inter.coupling.matrix{1,1}=eeqq2nqi(3.2e6,0.16,5/2,[0 0 0]);
inter.zeeman.eigs={[-5 -5 10]};
inter.zeeman.euler={[0 0 0]};

% Basis set
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spin system
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% MAS parameters and central-transition observable
parameters.spins={'27Al'};
parameters.rframes={{'27Al',3}};
rho_ct=sparse(3,4,1,6,6);
parameters.rho0=rho_ct(:);
parameters.coil=parameters.rho0;
parameters.ref_norm=norm(parameters.rho0)*norm(parameters.coil);
parameters.rf_op=operator(spin_system,'Lx','27Al');
parameters.rf_amp=0;

% Resolve the third-order central-transition evolution on both routes
[fp_nqi,sl_nqi]=compare_routes(spin_system,parameters,...
                               [3 5 7 9],[17 33 65 129]);
check_limit('27Al Q3 central transition',fp_nqi,sl_nqi,1e-4);

% Confirm that the third-order correction contributes
parameters.max_rank=1;
parameters.masframe='rotor';
parameters.orientation=[0 0 0];
third_gen=rotor_stack(spin_system,parameters,'labframe');
parameters.rframes={{'27Al',2}};
second_gen=rotor_stack(spin_system,parameters,'labframe');
q3_norm=norm(third_gen{1}-second_gen{1},'fro');
assert(q3_norm>100*eps(norm(third_gen{1},'fro')),...
       'The third-order quadrupolar correction was not resolved.');
fprintf('Third-order correction norm: %.9g rad/s\n',q3_norm);
fprintf('MAS_FPLANCK_SLICES_SUCCESS\n');

end

% Compare phase-space propagation with a midpoint rotor stack
function [fp_sig,sl_sig]=compare_routes(spin_system,parameters,...
                                         rotor_ranks,slice_counts)

% Propagate a single-crystal phase delta with the FP rotor generator
fp_sig=zeros(size(rotor_ranks));
for rank_idx=1:numel(rotor_ranks)
    parameters.max_rank=rotor_ranks(rank_idx);
    fp_sig(rank_idx)=singlerot(spin_system,@fp_signal,parameters,'labframe');
    fprintf('FP rank %d: %.9g%+.9gi\n',rotor_ranks(rank_idx),...
            real(fp_sig(rank_idx)),imag(fp_sig(rank_idx)));
end

% Traverse midpoint slices toward decreasing phase as in the FP generator
sl_sig=zeros(size(slice_counts));
for slice_idx=1:numel(slice_counts)
    slice_count=slice_counts(slice_idx);
    parameters.max_rank=(slice_count-1)/2;
    parameters.orientation=[0 0 -pi/slice_count];
    parameters.masframe='rotor';
    liouv_stack=rotor_stack(spin_system,parameters,'labframe');
    rho=parameters.rho0;
    for step_idx=1:slice_count
        rotor_idx=mod(1-step_idx,slice_count)+1;
        rho=expm(-1i*full(liouv_stack{rotor_idx}+...
                          parameters.rf_amp*parameters.rf_op)*...
                 (parameters.duration/slice_count))*rho;
    end
    sl_sig(slice_idx)=(parameters.coil'*rho)/parameters.ref_norm;
    fprintf('Slices %d: %.9g%+.9gi\n',slice_count,...
            real(sl_sig(slice_idx)),imag(sl_sig(slice_idx)));
end

end

% Add the same transverse RF operator at every FP rotor collocation point
function signal=fp_signal(~,parameters,generator,~,~)

% Add the RF operator at every rotor collocation point
generator=generator+parameters.rf_amp*...
          kron(speye(parameters.spc_dim),parameters.rf_op);

% Propagate the single-crystal phase state for one period
rho=expm(-1i*full(generator)*parameters.duration)*parameters.rho0;

% Detect and normalise the signal
signal=(parameters.coil'*rho)/parameters.ref_norm;

end

% Require both independent refinements to meet a normalised-signal tolerance
function check_limit(label,fp_sig,sl_sig,tolerance)

% Test the last refinement and the cross-route complex-signal difference
fp_step=abs(fp_sig(end)-fp_sig(end-1));
sl_step=abs(sl_sig(end)-sl_sig(end-1));
route_gap=abs(fp_sig(end)-sl_sig(end));
fprintf('%s: FP step %.6g, slice step %.6g, gap %.6g\n',...
        label,fp_step,sl_step,route_gap);
assert(all(isfinite([fp_step sl_step route_gap]))&&...
       fp_step<tolerance/2&&sl_step<tolerance/2&&route_gap<tolerance,...
       '%s MAS routes have not converged to the same signal.',label);

end


