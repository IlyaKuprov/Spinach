% Validate externally supplied reduced Lunghi two-proton fixtures.
%
%                    result=test_cce_lunghi_fixtures()
%
% Set SPINACH_CCE_LUNGHI_FIXTURES to inputs and CSVs, and
% SPINACH_CCE_LUNGHI_STATES to phase-preserving *_states.json files and
% *_diagnostic.mat source matrices exported by PyCCE. Source
% structures are not redistributed: the pinned source has no licence.
% These reduced source-input cases are not published-figure validation.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_lunghi_fixtures()

% Keep third-party inputs external until redistribution permission exists
result=new_test_result('kernel/cce_lunghi_fixtures','CCE Lunghi fixtures',...
                       'Reduced VOTPP and Cu(mnt)2 source-input gCCE reconstruction.');
folder=getenv('SPINACH_CCE_LUNGHI_FIXTURES'); states_folder=getenv('SPINACH_CCE_LUNGHI_STATES');
if isempty(folder)||isempty(states_folder)
    result.status='SKIP'; result.messages={'External unlicensed fixture directory not configured.'}; return;
end
errors=zeros(1,2); cases={'VOTPP','Cumnt'};
for case_id=1:2

    % Read explicit positions, anisotropic gyros and frequency tensors
    name=cases{case_id}; input=jsondecode(fileread(fullfile(folder,[name '_input.json'])));
    expected=readmatrix(fullfile(folder,[name '_gcce2.csv']));
    dims=2*input.center.s(:).'+1; central_dimension=prod(dims);
    if dims(1)==8, isotope='51V'; else, isotope='63Cu'; end
    sys.isotopes={isotope,'E','1H','1H'}; sys.magnet=0;
    sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
    inter.zeeman.matrix=repmat({zeros(3)},1,4); inter.coupling.matrix=cell(4);
    ss=create(sys,inter); field=input.simulator.magnetic_field(:).';
    gyros=cell(1,4); xyz=[input.center.xyz;input.bath.xyz];
    for j=1:2, gyros{j}=squeeze(input.center.gyro(j,:,:)); end
    for j=3:4, gyros{j}=input.bath.gyro(j-2)*eye(3); end

    % Match scalar nuclear and matrix central Zeeman conventions explicitly
    for j=1:4
        tensor=zeros(3);
        if j==2, tensor(3,:)=1000*field*gyros{j};
        else, tensor(3,:)=-1000*field*gyros{j}; end
        ss.inter.zeeman.matrix{j}=tensor;
    end
    d=input.center.D_kHz(1); e=input.center.D_kHz(2);
    ss.inter.coupling.matrix{1,1}=2*pi*1000*diag([-d/3+e -d/3-e 2*d/3]);
    ss.inter.coupling.matrix{1,2}=2*pi*1000*input.center.A_kHz.';

    % Construct source-convention point dipoles independently of the CCE builder
    for j=1:4
        for k=max(j+1,3):4
            r=xyz(j,:)-xyz(k,:); distance=norm(r);
            tensor=1.05457172*(eye(3)-3*(r.'*r)/distance^2)/distance^3;
            tensor=gyros{j}*tensor*gyros{k};
            ss.inter.coupling.matrix{j,k}=1000*tensor.';
        end
    end

    % Derive the logical central eigenstates and selective Hahn rotation
    p.central_spins=[1 2]; p.bath_spins=[3 4]; p.method='gcce';
    p.expansion='standard'; p.quantity='coherence'; p.order=2;
    p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
    p.graph=struct('mode','complete','score','coupling_norm'); p.meanfield='none';
    p.bath_state.mode='mixed_product'; p.times=input.compute.times_ms(:).'*1e-3;
    p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);
    p.sequence=[]; p.central.alpha=eye(central_dimension,1);
    p.central.beta=circshift(p.central.alpha,1);
    model=cce_prepare(ss,p); [hamiltonian,~]=cce_hamiltonian(model,[],struct());
    states=jsondecode(fileread(fullfile(states_folder,[name '_states.json'])));
    a=states.alpha_real+1i*states.alpha_imag; b=states.beta_real+1i*states.beta_imag;
    source=load(fullfile(states_folder,[name '_diagnostic.mat']),'H','Hc');
    central_error=norm(hamiltonian/(2*pi*1000)-source.Hc,'fro');
    matrix_scale=norm(source.Hc,'fro');
    result=test_close(result,[name ' central source Hamiltonian'],...
                      full(hamiltonian)/(2*pi*1000)/matrix_scale,source.Hc/matrix_scale,1e-10,0,...
                      'Explicit kHz to rad/s boundary and central Hilbert ordering.');
    residual=max(norm(hamiltonian*a-a*(a'*hamiltonian*a)),...
                 norm(hamiltonian*b-b*(b'*hamiltonian*b)))/norm(hamiltonian,'fro');
    result=test_true(result,[name ' exported eigenstates'],residual<100*eps,...
                     'Explicit source eigenvector phases define initial state and logical axes.');
    p.central.alpha=a; p.central.beta=b; p.central.rho=(a+b)*(a+b)'/2;
    p.sequence=struct('type','ideal','target',[1 2],'placement','fraction','at',0.5,...
                      'operator',pi*(a*b'+b*a')/2);
    p.normalized=true; actual=cce(ss,p);

    % Separate the stored source pulse from the unitary Spinach logical pulse
    [full_hamiltonian,~]=cce_hamiltonian(model,[3 4],struct());

    % Permute bath-first PyCCE ordering into the central-first native basis
    permutation=reshape(1:4*central_dimension,central_dimension,4).';
    permutation=permutation(:); reference_h=source.H(permutation,permutation);
    matrix_error=norm(full(full_hamiltonian)/(2*pi*1000)-reference_h,'fro');
    matrix_scale=norm(reference_h,'fro');
    result=test_close(result,[name ' complete source Hamiltonian'],...
                      full(full_hamiltonian)/(2*pi*1000)/matrix_scale,reference_h/matrix_scale,1e-10,0,...
                      'Independent PyCCE tensor compilation and explicit basis permutation.');
    fprintf('%s source matrices: absolute %.6g/%.6g kHz; relative %.6g/%.6g.\n',...
            name,central_error,matrix_error,central_error/norm(source.Hc,'fro'),...
            matrix_error/norm(reference_h,'fro'));

    % Propagate source and unitary pulses with dense matrix exponentials
    initial=kron((a+b)*(a+b)'/2,eye(4)/4);
    source_pulse=-1i*kron(a*b'+b*a',eye(4));
    unitary_pulse=expm(-1i*pi*kron(a*b'+b*a',eye(4))/2);
    source=zeros(numel(p.times),1); independent=source;
    for k=1:numel(p.times)
        free=expm(-1i*full(full_hamiltonian)*p.times(k)/2);
        evolution=free*source_pulse*free; density=evolution*initial*evolution';
        reduced=partial_trace(density,central_dimension);
        source(k)=2*b'*reduced*a;
        evolution=free*unitary_pulse*free; density=evolution*initial*evolution';
        reduced=partial_trace(density,central_dimension);
        independent(k)=2*a'*reduced*b;
    end
    reference=expected(:,2)+1i*expected(:,3); errors(case_id)=max(abs(source-reference));
    result=test_close(result,[name ' stored source pulse'],source,reference,1e-10,0,...
                      'Stored source uses -i sigma_x with zero spectator block.');
    result=test_close(result,[name ' unitary CCE pulse'],actual.value,independent,1e-10,0,...
                      'Spinach uses exp(-i pi sigma_x/2), retaining spectator levels.');
    result=test_true(result,[name ' distinct pulse semantics'],norm(source_pulse'*source_pulse-eye(4*central_dimension),'fro')>1,...
                     'Do not equate the source fixture pulse with a full unitary.');
end
disp(struct('VOTPP_error',errors(1),'Cumnt_error',errors(2)));

end

function reduced=partial_trace(density,dimension)

% Sum bath diagonal blocks in central-first Hilbert order
reduced=zeros(dimension);
for k=1:4, reduced=reduced+density(k:4:end,k:4:end); end

end


