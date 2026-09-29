% Compares explicitly sliced MAS evolution with Fokker-Planck Liouville MAS.
% Syntax:
%
%                     mas_fplanck_slices()
%
% No inputs. Prints convergence metrics; throws on a failed equivalence test.
% The cases are 13C CSA under finite RF and the central transition of a
% strongly quadrupolar 27Al spin with third-order rotating-frame terms.
%
% Calculation time: minutes
%
% talos@spindynamics.org

function mas_fplanck_slices()

% Set a single-spin CSA with RF and a phase-sensitive L+ to Lz transfer
sys.magnet=9.4;
sys.isotopes={'13C'};
sys.parallel={'processes',1};
sys.parprops={};
inter.zeeman.eigs={[-120 -25 145]};
inter.zeeman.euler={[0.4 0.7 0.2]};
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
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

% Refine the rotor rank and the midpoint slice count independently
[fp_csa,sl_csa]=compare_routes(spin_system,parameters,...
                               [4 6 8 10],[33 65 129]);
check_limit('13C CSA with RF',fp_csa,sl_csa,1e-3);

% Use the Smelko 27Al quadrupolar parameters with a selective central coherence
clear sys inter bas
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'27Al'};
sys.parallel={'processes',1};
sys.parprops={};
inter.coupling.matrix{1,1}=eeqq2nqi(3.2e6,0.16,5/2,[0 0 0]);
inter.zeeman.eigs={[-5 -5 10]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
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

% Check that third-order correction is exercised, not merely requested
parameters.max_rank=1;
parameters.masframe='rotor';
parameters.orientation=[0 0 0];
L3=rotor_stack(spin_system,parameters,'labframe');
parameters.rframes={{'27Al',2}};
L2=rotor_stack(spin_system,parameters,'labframe');
q3_norm=norm(L3{1}-L2{1},'fro');
assert(q3_norm>100*eps(norm(L3{1},'fro')),...
       'The third-order quadrupolar correction was not resolved.');
fprintf('Third-order correction norm: %.9g rad/s\n',q3_norm);
fprintf('MAS_FPLANCK_SLICES_SUCCESS\n');

end

% Compare phase-space propagation with a midpoint rotor stack
function [fp_sig,sl_sig]=compare_routes(spin_system,parameters,...
                                         ranks,counts)

% Propagate a single-crystal phase delta with the FP rotor generator
fp_sig=zeros(size(ranks));
for n=1:numel(ranks)
    parameters.max_rank=ranks(n);
    fp_sig(n)=singlerot(spin_system,@fp_signal,parameters,'labframe');
    fprintf('FP rank %d: %.9g%+.9gi\n',ranks(n),...
            real(fp_sig(n)),imag(fp_sig(n)));
end

% Traverse midpoint slices toward decreasing phase as in the FP generator
sl_sig=zeros(size(counts));
for n=1:numel(counts)
    count=counts(n);
    parameters.max_rank=(count-1)/2;
    parameters.orientation=[0 0 -pi/count];
    parameters.masframe='rotor';
    L=rotor_stack(spin_system,parameters,'labframe');
    rho=parameters.rho0;
    for k=1:count
        idx=mod(1-k,count)+1;
        rho=expm(-1i*full(L{idx}+parameters.rf_amp*parameters.rf_op)*...
                 (parameters.duration/count))*rho;
    end
    sl_sig(n)=(parameters.coil'*rho)/parameters.ref_norm;
    fprintf('Slices %d: %.9g%+.9gi\n',count,...
            real(sl_sig(n)),imag(sl_sig(n)));
end

end

% Add the same transverse RF operator at every FP rotor collocation point
function signal=fp_signal(~,parameters,G,~,~)

% Integrate over one rotor period and contract with the physical coil
G=G+parameters.rf_amp*kron(speye(parameters.spc_dim),parameters.rf_op);
rho=expm(-1i*full(G)*parameters.duration)*parameters.rho0;
signal=(parameters.coil'*rho)/parameters.ref_norm;

end

% Require both independent refinements to meet a normalised-signal target
function check_limit(label,fp_sig,sl_sig,target)

% Test the last refinement and the cross-route complex-signal difference
fp_step=abs(fp_sig(end)-fp_sig(end-1));
sl_step=abs(sl_sig(end)-sl_sig(end-1));
route_gap=abs(fp_sig(end)-sl_sig(end));
fprintf('%s: FP step %.6g, slice step %.6g, gap %.6g\n',...
        label,fp_step,sl_step,route_gap);
assert(all(isfinite([fp_step sl_step route_gap]))&&...
       fp_step<target/2&&sl_step<target/2&&route_gap<target,...
       '%s MAS routes have not converged to the same signal.',label);

end


