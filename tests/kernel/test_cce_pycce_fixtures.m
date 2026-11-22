% Test: kernel/cce_pycce_fixtures | PyCCE fixture parity
% Reproduce pinned PyCCE raw singleton/pair fixtures and reductions. Syntax:
%
%                    result=test_cce_pycce_fixtures()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% Data: MICCoMpy/PyCCE d74d3cb683bf3caa07231654314e554f440833b4.
% test_cce_pycce_data.json retains Stage-A inputs and complex raw traces.
% kHz, ms, Gauss, and rad/(ms Gauss) are converted explicitly to rad/s and s.
% Alpha/beta are source eigenvectors including their relative gauge. PyCCE
% observes the pulse-rotated observed basis and reverses Hahn manifold products
% in its projected engine. The projected source signal already matches the
% fixed receiver; only the generalised source signal is conjugated here.
% Logical-pi PyCCE projections annihilate spectator levels; our physical
% unitary preserves them. The discrepancy is separately measured, not hidden.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_pycce_fixtures()

% Load compact, source-pinned input and raw-output records
result=new_test_result('kernel/cce_pycce_fixtures','PyCCE fixture parity',...
                       'Three fixtures, raw subsets and complex final reductions.');
fixtures=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),...
                                      'test_cce_pycce_data.json')));
raw_errors=zeros(1,numel(fixtures));
final_errors=raw_errors;
bounds=raw_errors;
for f=1:numel(fixtures)
    fixture=fixtures(f);
    nb=numel(fixture.bath.N);

    % Initialise Spinach isotope metadata then install explicit canonical tensors
    sys.isotopes=[{'14N'} reshape(fixture.bath.N,1,[])];
    sys.magnet=0;
    sys.output='hush';
    sys.disable={'hygiene'};
    sys.parallel={'local',1};
    inter.zeeman.matrix=repmat({zeros(3)},1,nb+1);
    inter.coupling.matrix=cell(nb+1);
    spin_system=create(sys,inter);
    field=fixture.simulator.magnetic_field(:).';

    % External PyCCE central labels map to the observed-spin Hamiltonian
    spin_system.inter.zeeman.matrix{1}(3,:)=-fixture.central.gyro*field*1000;
    D=fixture.central.D;
    E=fixture.central.E;
    spin_system=cce_test_write(spin_system,[1 1],2*pi*1000*diag([-D/3+E -D/3-E 2*D/3]));
    gyros=zeros(1,nb);
    for j=1:nb
        type=fixture.bath.types.(matlab.lang.makeValidName(fixture.bath.N{j}));
        gyros(j)=type.gyro;
        spin_system.inter.zeeman.matrix{j+1}(3,:)=-type.gyro*field*1000;
        spin_system=cce_test_write(spin_system,[1 j+1],2*pi*1000*squeeze(fixture.bath.A(j,:,:)).');

        % Conventional reference CCE explicitly retains only observed Sz coupling
        if strcmp(fixture.compute.method,'cce')
            tensor_value=cce_test_read(spin_system,[1 j+1]);
            tensor_value(:,1:2)=0;
            spin_system=cce_test_write(spin_system,[1 j+1],tensor_value);
        end
        spin_system=cce_test_write(spin_system,[j+1 j+1],2*pi*1000*squeeze(fixture.bath.Q(j,:,:)).');
    end
    for j=1:nb
        for k=j+1:nb
            r=fixture.bath.xyz(j,:)-fixture.bath.xyz(k,:);
            distance=norm(r);
            tensor=1000*1.05457172*gyros(j)*gyros(k)*(eye(3)-3*(r.'*r)/distance^2)/distance^3;
            spin_system=cce_test_write(spin_system,[j+1 k+1],tensor.');
        end
    end

    % Resolve the exact observed-state and selective logical-pulse convention
    p.observed=1;
    p.cce.bath_spins=2:nb+1;
    p.cce.approximation='generalized';
    if strcmp(fixture.compute.method,'cce')
        p.cce.approximation='projected';
    end
    convention=-1;
    if strcmp(fixture.compute.method,'cce')
        convention=1;
    end
    p.cce.order=2;
    p.orientation=[0 0 0];
    p.cce.interaction_cutoff=Inf;
    p.cce.graph=struct('mode','complete','score','coupling_norm');
    p.cce.meanfield='none';
    p.cce.bath_state=struct('mode','mixed_product');
    p.times=fixture.times_ms(:).'/1000;
    p.cce.numerics=struct('atol',1e-5,'rtol',1e-5,'max_memory_bytes',2^29);
    a=fixture.alpha_real+1i*fixture.alpha_imag;
    b=fixture.beta_real+1i*fixture.beta_imag;
    p.rho0=(a+b)*(a+b)'/2;
    p.transition=struct('alpha',a,'beta',b);
    pulse_sequence=@hahn_echo;
    p.normalized=true;
    p.transition.spins=p.observed;
    [value,actual]=clexpn(spin_system,pulse_sequence,p,'labframe');

    % Compare the raw traces, not only the final absolute echo envelope
    for j=1:numel(fixture.raw)
        expected=fixture.raw(j);
        sites=sort(expected.sites(:).');
        index=find(arrayfun(@(r)isequal(r.sites,sites),actual.clusters),1)+1;
        if isempty(sites)
            index=1;
        end
        record=actual.raw{1}(index);
        reference=expected.real+1i*convention*expected.imag;
        discrepancy=max(abs(record.value-reference));
        raw_errors(f)=max(raw_errors(f),discrepancy);
        bound=max(record.abs_error_bound);
        bounds(f)=max(bounds(f),bound);
        result=test_true(result,'PyCCE raw subset',discrepancy<=2*bound,...
                         'Two eigensolvers plus exact decimal input serialization bounds.');
    end

    % Preserve external signal-column labels and the documented phase convention
    reference=fixture.coherence_real+1i*convention*fixture.coherence_imag;
    final_errors(f)=max(abs(value-reference));
    result=test_true(result,'PyCCE final signal',...
                     final_errors(f)<=2*max(actual.reduction.abs_error_bound),...
                     'Complex amplitude and phase are compared without clipping.');

    % Isolate the subset algebra using only the externally supplied raw traces
    records=actual.raw{1};
    records(1).value(:)=1;
    for j=1:numel(fixture.raw)
        expected=fixture.raw(j);
        sites=sort(expected.sites(:).');
        index=find(arrayfun(@(r)isequal(r.sites,sites),actual.clusters),1)+1;
        if isempty(sites)
            index=1;
        end
        records(index).value=expected.real+1i*convention*expected.imag;
    end
    for j=1:numel(records)
        records(j).abs_error_bound(:)=0;
    end
    reduction=struct('method',fixture.compute.method,'order',p.cce.order,...
                     'numerics',p.cce.numerics);
    reduction.quantity='signal';
    reduction.times=p.times;
    reduced=cce_reduce(records,actual.clusters,reduction);
    operations=numel(records)^2*numel(p.times);
    rounding=operations*eps/(1-operations*eps);
    result=test_close(result,'PyCCE reduction only',reduced.value,reference,...
                      rounding,rounding,'Same raw traces, independently compiled subset powers.');
end
metrics.raw_errors=raw_errors;
metrics.final_errors=final_errors;
metrics.maximum_raw_bounds=bounds;

% Report the measured metrics
disp(metrics);

end


