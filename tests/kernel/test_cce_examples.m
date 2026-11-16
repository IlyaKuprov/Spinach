% Smoke and independent-reference checks for the CCE application examples.
% Syntax: [result,reference_capture]=test_cce_examples()
% The optional capture retains independent physical arrays for regression.
%
% ilya.kuprov@weizmann.ac.il

function [result,reference_capture]=test_cce_examples()

% Initialise the standard test record and preserve figure visibility
result=new_test_result('kernel/cce_examples','CCE application examples',...
                       'Runnable categories and independently evaluated small limits.');
reference_capture=struct();
visibility=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off');

% Enumerate examples explicitly so new unrelated files cannot enter the test
names={'cce_single_nucleus','cce_nv_diamond','cce_sip_donor',...
       'cce_sic_divacancy','gcce_sic_clock','cce_p1_electron_bath',...
       'pcce_electron_layer','gcce_votpp','cce_erbium_cawo4',...
       'gcce_highspin_siv','cce_hbn_quadrupoles','cce_polarized_bath',...
       'gcce_two_centres','gcce_finite_drive','cce_deer_bandwidth',...
       'cce_noise_spectrum','mecce_dissipative_bath',...
       'cce_methyl_molecular','cce_molecular_powder'};
for n=1:numel(names)

    % Execute each public example and retain the standard harness fields only
    started=tic;
    if ismember(names{n},{'cce_noise_spectrum','cce_single_nucleus','cce_nv_diamond',...
                         'cce_p1_electron_bath','pcce_electron_layer',...
                         'cce_erbium_cawo4','cce_hbn_quadrupoles',...
                         'cce_methyl_molecular','cce_sic_divacancy','cce_sip_donor'})
        [answer,diagnostics]=feval(names{n});
        switch names{n}
            case 'cce_noise_spectrum'
                [reference_error,snapshot]=reference_cce_noise_spectrum(answer);
            case 'cce_single_nucleus'
                [reference_error,snapshot]=reference_cce_single_nucleus(answer);
            case 'cce_nv_diamond'
                [reference_error,snapshot]=reference_cce_nv_diamond(answer);
            case 'cce_p1_electron_bath'
                [reference_error,snapshot]=reference_cce_p1_electron_bath(answer);
            case 'pcce_electron_layer'
                [reference_error,snapshot]=reference_pcce_electron_layer(answer);
                result=test_true(result,'Partition full-order identity',...
                    snapshot.full_order_error<1e-9,...
                    'The retained pair agrees with the unpartitioned full bath.');
            case 'cce_erbium_cawo4'
                [reference_error,snapshot]=reference_cce_erbium_cawo4(answer);
            case 'cce_hbn_quadrupoles'
                [reference_error,snapshot]=reference_cce_hbn_quadrupoles(answer);
            case 'cce_methyl_molecular'
                [reference_error,snapshot]=reference_cce_methyl_molecular(answer);
            case 'cce_sic_divacancy'
                [reference_error,snapshot]=reference_cce_sic_divacancy(answer);
            case 'cce_sip_donor'
                [reference_error,snapshot]=reference_cce_sip_donor(answer);
        end
        actual=struct('value',answer,'reference_error',reference_error,...
                      'valid',diagnostics.reduction.valid,'density',[]);
        reference_capture.(names{n})=snapshot;
        reference_capture.(names{n}).diagnostics=diagnostics;
        result=test_true(result,[names{n} ' projected error bounds'],...
                         all(diagnostics.reduction.abs_error_bound>0),...
                         'Projected spectral error bounds must reach the reducer.');
    elseif ismember(names{n},{'gcce_sic_clock','gcce_highspin_siv','gcce_votpp',...
                             'gcce_two_centres','gcce_finite_drive','mecce_dissipative_bath'})
        [answer,diagnostics]=feval(names{n});
        actual=diagnostics.reduction; actual.times=diagnostics.parameters.times(:);

        % Check the observed density without repairing its populations
        dimension=sqrt(size(answer,2));
        actual.density=reshape(answer.',dimension,dimension,[]);
        actual.diagnostics.trace_error=zeros(size(answer,1),1);
        actual.diagnostics.minimum_eigenvalue=zeros(size(answer,1),1);
        for sample=1:size(answer,1)
            rho=actual.density(:,:,sample);
            actual.diagnostics.trace_error(sample)=abs(trace(rho)-1);
            actual.diagnostics.minimum_eigenvalue(sample)=min(eig((rho+rho')/2));
        end
        switch names{n}
            case 'gcce_sic_clock'
                [actual.reference_error,reference_capture.gcce_sic_clock]=gcce_sic_clock_reference(answer);
            case 'gcce_highspin_siv'
                [actual.reference_error,reference_capture.gcce_highspin_siv]=gcce_highspin_siv_reference(answer);
            case 'gcce_votpp'
                [actual.reference_error,reference_capture.gcce_votpp]=gcce_votpp_reference(answer);
            case 'gcce_finite_drive'
                [actual.reference_error,reference_capture.gcce_finite_drive]=gcce_finite_drive_reference(answer,diagnostics);
            case 'gcce_two_centres'
                [actual.reference_error,reference_capture.gcce_two_centres]=...
                    gcce_two_centres_reference(answer,diagnostics);
            case 'mecce_dissipative_bath'
                [actual.reference_error,reference_capture.mecce_dissipative_bath]=...
                    mecce_dissipative_bath_reference(answer,diagnostics);
                result=test_true(result,'Liouville arithmetic bounds',...
                    all(diagnostics.reduction.abs_error_bound>0,'all'),...
                    'The stock Liouville action retains its arithmetic estimate.');
        end
        reference_capture.(names{n}).diagnostics=diagnostics;
    elseif strcmp(names{n},'cce_polarized_bath')
        [answer,diagnostics]=cce_polarized_bath();
        actual=diagnostics.reduction; actual.value=answer;
        [actual.sampling_error,reference_capture.cce_polarized_bath]=...
            cce_polarized_bath_reference(answer,diagnostics);
        actual.density=[];
        reference_capture.cce_polarized_bath.diagnostics=diagnostics;
        parts=diagnostics.reduction.samples; weights=diagnostics.reduction.weights;
        bound=zeros(size(answer)); magnitude=zeros(size(answer));
        for sample=1:numel(parts)
            bound=bound+weights(sample)*parts{sample}.abs_error_bound;
            magnitude=magnitude+weights(sample)*abs(parts{sample}.value);
        end
        operations=4*numel(parts); rounding=operations*eps/(1-operations*eps);
        result=test_close(result,'Sampled bath arithmetic bound',...
            diagnostics.reduction.abs_error_bound,bound+rounding*magnitude,0,0,...
            'Deterministic summation error is distinct from estimator uncertainty.');
    else
        actual=feval(names{n});
    end
    result=test_true(result,[names{n} ' valid'],all(actual.valid,'all'),...
                     'Every displayed time point has a valid expansion.');
    result=test_true(result,[names{n} ' finite'],all(isfinite(actual.value),'all'),...
                     'Every returned complex value is finite.');
    if isfield(actual,'reference_error')

        % Independent small-system targets use an absolute trace-error criterion
        result=test_close(result,[names{n} ' reference'],actual.reference_error,...
                          0,1e-7,0,'Analytic, dense expm, or pinned Python reference.');
    end
    if ~isempty(actual.density)

        % Full-order finite baths must preserve observed density physicality
        result=test_true(result,[names{n} ' trace'],...
                         max(actual.diagnostics.trace_error)<1e-8,...
                         'Reduced density matrices retain unit trace.');
        result=test_true(result,[names{n} ' positive'],...
                         min(actual.diagnostics.minimum_eigenvalue)>-1e-8,...
                         'No material negative population at full bath order.');
    end
    if strcmp(names{n},'cce_deer_bandwidth')

        % Excluding the pump resonance must remove its selected pulse
        parameters=actual.input_parameters;
        parameters.deer.detuning_center=2*pi*1e7;
        sys.magnet=0.01; sys.isotopes={'E','E','E'}; sys.output='hush';
        inter.zeeman.scalar={0,0,0};
        model=cce_prepare(create(sys,inter),parameters);
        schedule=cce_controls(model,[],[],[],parameters);
        result=test_true(result,'DEER example excludes pump resonance',...
                         isempty(schedule.selected_sites{2}),...
                         'The pump pulse must follow the documented frequency window.');
    end
    if strcmp(names{n},'cce_molecular_powder')

        % All public response arrays must describe the powder average
        result=test_close(result,'Powder response alias',actual.response,...
                          actual.value,0,0,...
                          'The response alias must equal the plotted powder signal.');
        result=test_close(result,'Powder final order response',...
                          actual.by_order(:,:,end),actual.value,1e-14,...
                          1e-14,'The final order must be orientation-averaged.');
        expected=zeros(size(actual.abs_error_bound));
        for k=1:numel(actual.orientation_weights)
            expected=expected+actual.orientation_weights(k)*...
                              actual.orientation_abs_error_bound(:,k);
        end
        result=test_close(result,'Powder averaged numerical bounds',...
                          actual.abs_error_bound,expected,0,0,...
                          'Positive quadrature weights also average absolute bounds.');
        result=test_close(result,'Powder diagnostic numerical bounds',...
                          actual.diagnostics.numerical_error,actual.abs_error_bound,...
                          0,0,...
                          'Diagnostics must be recomputed for the powder result.');
    end
    fprintf('Example %s completed in %.3f s.\n',names{n},toc(started));
    close all;
end

end


% Independent dense reference retained from gcce_sic_clock
function [reference_error,snapshot]=gcce_sic_clock_reference(answer)

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E3','13C','29Si'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,1}=diag([-1.0e4+2e3 -1.0e4-2e3 2e4]);
inter.coupling.matrix{1,2}=diag([2200 1800 3000]);
inter.coupling.matrix{1,3}=diag([900 1200 1700]);
inter.coupling.matrix{2,3}=diag([120 120 -240]);

% Compile physical metadata and the independently specified preparation
spin_system=create(sys,inter);
v=[1;0;1]/sqrt(2); observed=v*v';
times=linspace(0,200e-6,31);

% Independently construct the full dense Hamiltonian from the input Hz tensors
dimensions=spin_system.comp.mults; dimension=prod(dimensions);
operators=cell(numel(dimensions),3); h=zeros(dimension);
for j=1:numel(dimensions)

    % Local angular momentum matrices and their tensor-product embeddings
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
for j=1:numel(dimensions)
    for k=j:numel(dimensions)

        % Spinach ordered pair matrices transpose into Cartesian operators
        tensor=inter.coupling.matrix{j,k};
        if isempty(tensor), continue; end
        for a=1:3
            for b=1:3
                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
            end
        end
    end
end

% Propagate an explicit product density and trace over the whole finite bath
d_obs=size(observed,1); db=dimension/d_obs;
rho0=kron(observed,eye(db)/db); reference=zeros(size(answer));
for n=1:numel(times)

    % This independent calculation uses expm, not the CCE eigensolver
    u=expm(-1i*h*times(n)); rho=u*rho0*u'; reduced=zeros(d_obs);
    for k=1:db
        ids=k:db:dimension; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
reference_error=max(abs(answer-reference),[],'all');

% Preserve independent physical arrays for the stage-0 comparison
snapshot=struct('a',a,'b',b,'db',db,'d_obs',d_obs,...
                'dimension',dimension,'dimensions',dimensions,'h',h,...
                'local',{local},'m',m,'operators',{operators},...
                'plus',plus,'reduced',reduced,'reference',reference,...
                'rho',rho,'rho0',rho0,'u',u,'v',v,'observed',observed,...
                'reference_error',reference_error,'tensor',tensor,'ids',ids);

end


% Independent dense reference retained from gcce_highspin_siv
function [reference_error,snapshot]=gcce_highspin_siv_reference(answer)

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E4','13C','29Si'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,1}=diag([-1 -1 2])*1e4;
inter.coupling.matrix{1,2}=diag([1700 1300 2400]);
inter.coupling.matrix{1,3}=diag([800 900 1500]);
inter.coupling.matrix{2,3}=diag([100 100 -200]);

% Compile physical metadata and the independently specified preparation
spin_system=create(sys,inter);
v=[0;1;1;0]/sqrt(2); observed=v*v';
times=linspace(0,200e-6,31);

% Independently construct the full dense Hamiltonian from the input Hz tensors
dimensions=spin_system.comp.mults; dimension=prod(dimensions);
operators=cell(numel(dimensions),3); h=zeros(dimension);
for j=1:numel(dimensions)

    % Local angular momentum matrices and their tensor-product embeddings
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
for j=1:numel(dimensions)
    for k=j:numel(dimensions)

        % Spinach ordered pair matrices transpose into Cartesian operators
        tensor=inter.coupling.matrix{j,k};
        if isempty(tensor), continue; end
        for a=1:3
            for b=1:3
                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
            end
        end
    end
end

% Propagate an explicit product density and trace over the whole finite bath
d_obs=size(observed,1); db=dimension/d_obs;
rho0=kron(observed,eye(db)/db); reference=zeros(size(answer));
for n=1:numel(times)

    % This independent calculation uses expm, not the CCE eigensolver
    u=expm(-1i*h*times(n)); rho=u*rho0*u'; reduced=zeros(d_obs);
    for k=1:db
        ids=k:db:dimension; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
reference_error=max(abs(answer-reference),[],'all');

% Preserve independent physical arrays for the stage-0 comparison
snapshot=struct('a',a,'b',b,'db',db,'d_obs',d_obs,...
                'dimension',dimension,'dimensions',dimensions,'h',h,...
                'local',{local},'m',m,'operators',{operators},...
                'plus',plus,'reduced',reduced,'reference',reference,...
                'rho',rho,'rho0',rho0,'u',u,'v',v,'observed',observed,...
                'reference_error',reference_error,'tensor',tensor,'ids',ids);

end


% Independent dense reference retained from gcce_votpp
function [reference_error,snapshot]=gcce_votpp_reference(answer)

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E','51V','1H'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,2}=diag([15 15 50])*1e6;
inter.coupling.matrix{1,3}=diag([2000 3000 7000]);
inter.coupling.matrix{2,3}=diag([100 100 -200]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*1e5;

% Compile physical metadata and the independently specified preparation
spin_system=create(sys,inter);
v=zeros(16,1); v([1 9])=1/sqrt(2); observed=v*v';
times=linspace(0,200e-6,31);

% Independently construct the full dense Hamiltonian from the input Hz tensors
dimensions=spin_system.comp.mults; dimension=prod(dimensions);
operators=cell(numel(dimensions),3); h=zeros(dimension);
for j=1:numel(dimensions)

    % Local angular momentum matrices and their tensor-product embeddings
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
for j=1:numel(dimensions)
    for k=j:numel(dimensions)

        % Spinach ordered pair matrices transpose into Cartesian operators
        tensor=inter.coupling.matrix{j,k};
        if isempty(tensor), continue; end
        for a=1:3
            for b=1:3
                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
            end
        end
    end
end

% Propagate an explicit product density and trace over the whole finite bath
d_obs=size(observed,1); db=dimension/d_obs;
rho0=kron(observed,eye(db)/db); reference=zeros(size(answer));
for n=1:numel(times)

    % This independent calculation uses expm, not the CCE eigensolver
    u=expm(-1i*h*times(n)); rho=u*rho0*u'; reduced=zeros(d_obs);
    for k=1:db
        ids=k:db:dimension; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
reference_error=max(abs(answer-reference),[],'all');

% Preserve independent physical arrays for the stage-0 comparison
snapshot=struct('a',a,'b',b,'db',db,'d_obs',d_obs,...
                'dimension',dimension,'dimensions',dimensions,'h',h,...
                'local',{local},'m',m,'operators',{operators},...
                'plus',plus,'reduced',reduced,'reference',reference,...
                'rho',rho,'rho0',rho0,'u',u,'v',v,'observed',observed,...
                'reference_error',reference_error,'tensor',tensor,'ids',ids);

end

% Independent dense projected reference for cce_erbium_cawo4
function [reference_error,snapshot]=reference_cce_erbium_cawo4(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','183W','183W'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[8 0 2];[0 10 -1]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;3000 1200 2000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-900 400 1000].';
inter.coupling.matrix{2,3}=diag([80 80 -160]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

spin_system=create(sys,inter);
parameters.bath=2:3;
parameters.times=linspace(0,200e-6,41);

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(parameters.times),1));
for n=1:numel(parameters.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*parameters.times(n)/2); ub=expm(-1i*hb*parameters.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
reference_error=max(abs(answer-reference));

snapshot=struct();
snapshot.a=a; snapshot.b=b;
snapshot.bath=bath;
snapshot.dimensions=dimensions;
snapshot.dimension=dimension;
snapshot.operators=operators;
snapshot.q=q;
snapshot.m=m;
snapshot.plus=plus;
snapshot.local=local;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.id=id;
snapshot.tensor=tensor;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.reference=reference;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;

end

% Independent dense projected reference for cce_hbn_quadrupoles
function [reference_error,snapshot]=reference_cce_hbn_quadrupoles(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','11B','14N'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[1.45 0 0];[-0.725 1.2557 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;40000 0 20000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;20000 0 -10000].';
inter.coupling.matrix{2,3}=diag([500 500 -1000]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
inter.coupling.matrix{3,3}=diag([-1 -1 2])*5e4;

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

spin_system=create(sys,inter);
parameters.bath=2:3;
parameters.times=linspace(0,200e-6,41);

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(parameters.times),1));
for n=1:numel(parameters.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*parameters.times(n)/2); ub=expm(-1i*hb*parameters.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
reference_error=max(abs(answer-reference));

snapshot=struct();
snapshot.a=a; snapshot.b=b;
snapshot.bath=bath;
snapshot.dimensions=dimensions;
snapshot.dimension=dimension;
snapshot.operators=operators;
snapshot.q=q;
snapshot.m=m;
snapshot.plus=plus;
snapshot.local=local;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.id=id;
snapshot.tensor=tensor;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.reference=reference;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;

end

% Independent dense projected reference for cce_methyl_molecular
function [reference_error,snapshot]=reference_cce_methyl_molecular(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','1H','1H','1H'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[5 1 0];[5 -0.5 0.866];[5 -0.5 -0.866]};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;7000 0 5000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;5000 0 4000].';
inter.coupling.matrix{1,4}=[0 0 0;0 0 0;6000 0 3000].';
inter.coupling.matrix{2,3}=2000*eye(3);
inter.coupling.matrix{2,4}=2000*eye(3);
inter.coupling.matrix{3,4}=2000*eye(3);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

spin_system=create(sys,inter);
parameters.bath=2:4;
parameters.times=linspace(0,200e-6,41);

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(parameters.times),1));
for n=1:numel(parameters.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*parameters.times(n)/2); ub=expm(-1i*hb*parameters.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
reference_error=max(abs(answer-reference));

snapshot=struct();
snapshot.a=a; snapshot.b=b;
snapshot.bath=bath;
snapshot.dimensions=dimensions;
snapshot.dimension=dimension;
snapshot.operators=operators;
snapshot.q=q;
snapshot.m=m;
snapshot.plus=plus;
snapshot.local=local;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.id=id;
snapshot.tensor=tensor;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.reference=reference;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;

end

% Independent dense projected reference for cce_sic_divacancy
function [reference_error,snapshot]=reference_cce_sic_divacancy(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','13C','29Si'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[4 0 1];[0 5 2]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 -4000].';
inter.coupling.matrix{2,3}=diag([250 250 -500]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

spin_system=create(sys,inter);
parameters.bath=2:3;
parameters.times=linspace(0,200e-6,41);

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(parameters.times),1));
for n=1:numel(parameters.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*parameters.times(n)/2); ub=expm(-1i*hb*parameters.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
reference_error=max(abs(answer-reference));

snapshot=struct();
snapshot.a=a; snapshot.b=b;
snapshot.bath=bath;
snapshot.dimensions=dimensions;
snapshot.dimension=dimension;
snapshot.operators=operators;
snapshot.q=q;
snapshot.m=m;
snapshot.plus=plus;
snapshot.local=local;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.id=id;
snapshot.tensor=tensor;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.reference=reference;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;

end

% Independent dense projected reference for cce_sip_donor
function [reference_error,snapshot]=reference_cce_sip_donor(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','29Si','29Si'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[8 4 0];[12 2 4]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
a=25; k0=0.85*2*pi/5.43;
for n=2:3

    % Hydrogenic envelope and six-valley interference, illustrative scale
    r=inter.coordinates{n};
    contact=2e5*exp(-2*norm(r)/a)*(sum(cos(k0*r))/3)^2;

    % Retain observed Sz only for the explicitly projected donor model
    inter.coupling.matrix{1,n}=diag([0 0 contact]);
end
inter.coupling.matrix{2,3}=diag([200 200 -400]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

spin_system=create(sys,inter);
parameters.bath=2:3;
parameters.times=linspace(0,200e-6,41);

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(parameters.times),1));
for n=1:numel(parameters.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*parameters.times(n)/2); ub=expm(-1i*hb*parameters.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
reference_error=max(abs(answer-reference));

snapshot=struct();
snapshot.a=a; snapshot.b=b;
snapshot.bath=bath;
snapshot.dimensions=dimensions;
snapshot.dimension=dimension;
snapshot.operators=operators;
snapshot.q=q;
snapshot.m=m;
snapshot.plus=plus;
snapshot.local=local;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.id=id;
snapshot.tensor=tensor;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.reference=reference;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;
snapshot.r=r; snapshot.contact=contact; snapshot.k0=k0;

end

% Independent reference retained from cce_single_nucleus
function [reference_error,snapshot]=reference_cce_single_nucleus(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E3','13C'};
sys.output='hush'; sys.parallel={'local',1};

% Nuclear Zeeman frequency and projected hyperfine components, rad/s
omega=-spin('13C')*sys.magnet; apar=0.3*abs(omega);
aper=0.6*abs(omega);

% Standard Spinach tensors: electron g, nuclear ppm, and coupling Hz
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=[0 0 aper;0 0 0;0 0 apar]/(2*pi);

% Compile metadata without constructing a full-bath basis

create(sys,inter);

% Projected signal; static explicit operators include their angle
parameters.observed_spins=1; parameters.bath_spins=2;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='signal'; parameters.order=1;
parameters.transition.alpha=[0;1;0]; parameters.transition.beta=[1;0;0];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,8*pi/abs(omega),101);
parameters.sequence=struct('type','ideal','target',1,...
                           'operator',pi*[0 1 0;1 0 0;0 0 0]/2,...
                           'angle',1,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^28);

% Run CCE and calculate an independent analytic echo
result=struct('value',answer,'times',parameters.times(:),'valid',true(size(answer)));
omega1=sqrt((omega+apar)^2+aper^2);
reference=1-2*(aper/omega1)^2*sin(omega*parameters.times(:)/4).^2.*...
          sin(omega1*parameters.times(:)/4).^2;
result.reference=reference; result.reference_error=max(abs(result.value-reference));
fprintf('Single-nucleus analytic Hahn maximum error: %.3g\n',result.reference_error);

reference_error=result.reference_error;
snapshot=struct();
snapshot.apar=apar;
snapshot.aper=aper;
snapshot.omega=omega;
snapshot.omega1=omega1;
snapshot.reference=reference;
snapshot.reference_error=reference_error;

end

% Independent reference retained from cce_nv_diamond
function [reference_error,snapshot]=reference_cce_nv_diamond(answer)

% Standard Spinach system with source gyromagnetic ratios mapped to g/ppm
sys.magnet=500*1e-4; sys.isotopes={'E3','13C','13C','13C','13C'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=cell(1,5); inter.coupling.matrix=cell(5);
gfree=2.00231930436256;
inter.zeeman.matrix{1}=gfree*(-17608.59705*1e7/spin('E3'))*eye(3);
for n=2:5

    % Source rad/(ms Gauss) to Spinach nuclear chemical shift in ppm
    inter.zeeman.matrix{n}=1e6*(6.7282853242747445*1e7/spin('13C')-1)*eye(3);
end

% Source zero-field splitting and explicit hyperfine tensors, converted kHz to Hz
inter.coupling.matrix{1,1}=1000*diag([-960000 -960000 1920000]);
inter.coupling.matrix{1,2}=1000*[62.09198322383363 42.918773452863142 86.078663610517651;42.918773452863142 -52.625936690744069 28.63925630406008;86.078663610517651 28.63925630406008 -9.4660465330895676];
inter.coupling.matrix{1,3}=1000*[46.846946007935479 41.7131493536275 20.798152899007547;41.7131493536275 -15.576682609001635 10.413661409461005;20.798152899007547 10.413661409461005 -31.270263398933857];
inter.coupling.matrix{1,4}=1000*[7.7040985952493219 -41.664463528969883 -27.802268186621028;-41.664463528969883 7.7040985952493219 27.802268186621028;-27.802268186621028 27.802268186621028 -15.408197190498676];
inter.coupling.matrix{1,5}=1000*[14.689630582202209 -18.600579331667948 -27.874817765944968;-18.600579331667948 -13.146092067661488 13.956956441013119;-27.874817765944968 13.956956441013119 -1.5435385145407345];

% Declare the reference's observed-Sz projected approximation before propagation
for n=2:5, inter.coupling.matrix{1,n}(:,1:2)=0; end

% Explicit bath coordinates and point-dipole bath couplings in Hz
xyz=[5.3499999999999996 1.78 3.5699999999999998;7.1299999999999999 3.5699999999999998 1.78;-5.3499999999999996 5.3499999999999996 3.5699999999999998;-7.1299999999999999 3.5699999999999998 5.3499999999999996];
for j=1:4
    for k=j+1:4

        % Independent source-unit conversion; displacement is in Angstrom
        r=(xyz(k,:)-xyz(j,:))*1e-10; u=r/norm(r);
        angular=1e-7*1.054571817e-34*(6.7282853242747445*1e7)^2/norm(r)^3;
        inter.coupling.matrix{j+1,k+1}=angular*(eye(3)-3*(u.'*u))/(2*pi);
    end
end


create(sys,inter);

% The source selected two-level transition is mS=-1/+1; a selective static operator includes its angle
parameters.observed_spins=1; parameters.bath_spins=2:5;
parameters.method='cce'; parameters.expansion='standard';
parameters.quantity='signal'; parameters.order=2;
parameters.transition.alpha=[0 0 1].';
parameters.transition.beta=[1 0 0].';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=1e-3*[0 0.012500000000000001 0.025000000000000001 0.037500000000000006 0.050000000000000003 0.0625 0.075000000000000011 0.087500000000000008 0.10000000000000001 0.1125 0.125 0.13750000000000001 0.15000000000000002 0.16250000000000001 0.17500000000000002 0.1875 0.20000000000000001 0.21250000000000002 0.22500000000000001 0.23750000000000002 0.25 0.26250000000000001 0.27500000000000002 0.28750000000000003 0.30000000000000004 0.3125 0.32500000000000001 0.33750000000000002 0.35000000000000003 0.36250000000000004 0.375 0.38750000000000001 0.40000000000000002 0.41250000000000003 0.42500000000000004 0.4375 0.45000000000000001 0.46250000000000002 0.47500000000000003 0.48750000000000004 0.5];
a=parameters.transition.alpha; b=parameters.transition.beta;
parameters.sequence=struct('type','ideal','target',1,'operator',...
                           pi*(a*b'+b*a')/2,'angle',1,...
                           'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-3,'rtol',1e-3,'max_memory_bytes',2^28);

% PyCCE reports the opposite off-diagonal element: conjugate its complex trace
result=struct('value',answer,'times',parameters.times(:),'valid',true(size(answer)));
reference=[1.0000000000000073 0.80412187352127196 0.97076472650584211 0.97345787001549855 0.92208749592281747 0.89487943408951576 0.97304683857265373 0.93078025774994988 0.86125671741742393 0.99661776302799532 0.76355437703758766 0.99904105154528788 0.9882134199268271 0.95383133371749873 0.86766134502112169 0.98965017367005748 0.89372993555346014 0.8732162021890324 0.98602797149751487 0.82356170821295749 0.99419409116745483 0.93938769167121006 0.9940891681376085 0.86777114512359543 0.9860547697233486 0.84220648442521762 0.93028570851733106 0.98826030630667416 0.82912496114994416 0.97421567083728422 0.95985232702599865 0.98902742105821617 0.86305059097592329 0.99445386679999404 0.79888714808310557 0.97200285714316348 0.96135900760617921 0.85148111754718481 0.96346493820569012 0.97793634578837296 0.9615196189824301].'+...
          1i*[3.4800759960904551e-33 -2.3669791480982818e-06 3.9053578383864199e-08 4.1251546477359897e-06 1.3322998842642245e-06 -1.9800296667105094e-05 2.7234071561477375e-06 -4.4544694809084375e-06 -4.2527858745940881e-06 -4.0487554818223102e-08 3.5977740510985286e-05 1.2251660190342334e-07 6.6845336712676628e-07 8.2300776394647603e-06 -1.6303824478936545e-05 9.7554193182357096e-07 -3.5152718524314725e-05 -4.581728372659567e-05 3.2348614343623995e-07 8.2517282974561609e-06 -2.5910341445473197e-06 2.4206767247835291e-05 3.6414673240542315e-06 1.0680014272253585e-05 5.5858781900956123e-06 4.8774427328959467e-05 1.3601894267783467e-05 -5.3571546260425495e-06 -1.8314667050526693e-05 -2.0325127606064874e-05 1.7291666493708969e-05 -2.6488761044315708e-06 -0.00015685884541566557 4.8639994138961556e-06 0.00013334140466861785 -9.8320747276684197e-06 3.1600752332004179e-05 0.00014601325410634766 9.4604540883105856e-06 1.1986705000717478e-05 -3.3514391192247085e-05].';
result.reference_error=max(abs(result.value-conj(reference)));
fprintf('NV pinned PyCCE CCE2 maximum complex error: %.3g.\n',result.reference_error);
fprintf('Seo ensemble T2=0.796 ms,n=2.220 is not a four-nucleus target.\n');

reference_error=result.reference_error;
snapshot=struct();
snapshot.a=a;
snapshot.angular=angular;
snapshot.b=b;
snapshot.gfree=gfree;
snapshot.r=r;
snapshot.reference=reference;
snapshot.u=u;
snapshot.xyz=xyz;
snapshot.reference_error=reference_error;

end

% Independent reference retained from cce_p1_electron_bath
function [reference_error,snapshot]=reference_cce_p1_electron_bath(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','E','14N'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[35 0 0];[35 0 1]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=diag([2000 2000 -4000]);
inter.coupling.matrix{2,3}=diag([81 81 114])*1e6;
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;400 0 200].';

% Declare observed-Sz projection for this projected electron-bath model
for n=2:numel(sys.isotopes)
    if ~isempty(inter.coupling.matrix{1,n}), inter.coupling.matrix{1,n}(:,1:2)=0; end
end

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Projected observed-spin Hahn signal and a complete small-bath graph
parameters.observed_spins=1; parameters.bath_spins=2:3;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='signal'; parameters.order=2;
parameters.transition.alpha=[1;0]; parameters.transition.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,41);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-3,'rtol',1e-3,'max_memory_bytes',2^28);

% Keep the electron and its nitrogen in the same physical group
parameters.expansion='partition'; parameters.groups={[2 3]};
parameters.partition=struct('mode','user','max_dimension',6);
parameters.order=1;

% Run the cluster expansion and report its physical input boundary
result=struct('value',answer,'times',parameters.times(:),'valid',true(size(answer)));
fprintf('cce_p1_electron_bath: %d/%d valid values; finite-cluster demonstration.\n',...
        nnz(result.valid),numel(result.valid));
fprintf('No P1 density/T2 claim: one explicit hyperfine-bound group.\n');

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath_spins; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(result.times),1));
for n=1:numel(result.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*result.times(n)/2); ub=expm(-1i*hb*result.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
result.reference_error=max(abs(result.value-reference));
fprintf('Independent dense Hahn maximum complex error: %.3g.\n',result.reference_error);

reference_error=result.reference_error;
snapshot=struct();
snapshot.a=a;
snapshot.b=b;
snapshot.bath=bath;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.dimension=dimension;
snapshot.dimensions=dimensions;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.id=id;
snapshot.local=local;
snapshot.m=m;
snapshot.operators=operators;
snapshot.plus=plus;
snapshot.reference=reference;
snapshot.tensor=tensor;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;

end

% Independent reference retained from pcce_electron_layer
function [reference_error,snapshot]=reference_pcce_electron_layer(answer)

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','E','E','E'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[20 0 0];[21 0 0];[0 30 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-12000 0 6000].';
inter.coupling.matrix{1,4}=[0 0 0;0 0 0;4000 0 3000].';
inter.coupling.matrix{2,3}=diag([2e5 2e5 -4e5]);
inter.coupling.matrix{3,4}=diag([300 300 -600]);

% Declare observed-Sz projection for this projected electron-bath model
for n=2:numel(sys.isotopes)
    if ~isempty(inter.coupling.matrix{1,n}), inter.coupling.matrix{1,n}(:,1:2)=0; end
end

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Projected observed-spin Hahn signal and a complete small-bath graph
parameters.observed_spins=1; parameters.bath_spins=2:4;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='signal'; parameters.order=3;
parameters.transition.alpha=[1;0]; parameters.transition.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,41);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-5,'rtol',1e-5,'max_memory_bytes',2^28);

% Freeze the strongly coupled pair as one partition
parameters.expansion='partition'; parameters.groups={[2 3],4};
parameters.partition=struct('mode','user','max_dimension',4);
parameters.order=2;

% Run the cluster expansion and report its physical input boundary
result=struct('value',answer,'times',parameters.times(:),'valid',true(size(answer)));
fprintf('pcce_electron_layer: %d/%d valid values; finite-cluster demonstration.\n',...
        nnz(result.valid),numel(result.valid));

% Compare with a separate unpartitioned full-order expansion
full_parameters=rmfield(parameters,{'partition','groups'});
full_parameters.expansion='standard'; full_parameters.order=3;
native_parameters=struct('observed',1,'times',parameters.times,...
    'transition',struct('spins',1,'ms',[1/2 -1/2]),...
    'cce',struct('bath_spins',2:4,'order',3,'approximation','projected',...
    'graph',parameters.graph,'numerics',parameters.numerics));
[full_result.value,full_diagnostics]=clexpn(spin_system,@hahn_echo,native_parameters,'labframe');
result.reference_error=max(abs(result.value-full_result.value));
fprintf('Partition/full-order maximum difference: %.3g (same Hamiltonian).\n',...
        result.reference_error);

% Independently assemble projected bath Hamiltonians from input Hz tensors
bath=parameters.bath_spins; dimensions=spin_system.comp.mults(bath);
dimension=prod(dimensions); operators=cell(numel(bath),3);
for j=1:numel(bath)

    % Angular momentum matrices in descending magnetic-quantum-number order
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3

        % Embed a physical bath operator with dense Kronecker products
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
bath_h=zeros(dimension); difference=zeros(dimension);
for j=1:numel(bath)

    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
    id=bath(j);
    if sys.isotopes{id}(1)~='E'
        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
    end
    for a=1:3

        % Projected hyperfine difference and on-site quadratic terms
        if ~isempty(inter.coupling.matrix{1,id})
            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
        end
        for b=1:3
            if ~isempty(inter.coupling.matrix{id,id})
                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
            end
        end
    end
    for k=j+1:numel(bath)
        for a=1:3
            for b=1:3

                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                tensor=inter.coupling.matrix{id,bath(k)};
                if ~isempty(tensor)
                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                end
            end
        end
    end
end
ha=bath_h+difference/2; hb=bath_h-difference/2;
reference=complex(zeros(numel(result.times),1));
for n=1:numel(result.times)

    % Independent full-bath Hahn propagators for a maximally mixed state
    ua=expm(-1i*ha*result.times(n)/2); ub=expm(-1i*hb*result.times(n)/2);
    reference(n)=trace(ub*ua*(ua*ub)')/dimension;
end
result.reference_error=max(abs(result.value-reference));
fprintf('Independent dense Hahn maximum complex error: %.3g.\n',result.reference_error);

reference_error=result.reference_error;
snapshot=struct();
snapshot.a=a;
snapshot.b=b;
snapshot.bath=bath;
snapshot.bath_h=bath_h;
snapshot.difference=difference;
snapshot.dimension=dimension;
snapshot.dimensions=dimensions;
snapshot.full_parameters=full_parameters;
snapshot.ha=ha;
snapshot.hb=hb;
snapshot.id=id;
snapshot.local=local;
snapshot.m=m;
snapshot.operators=operators;
snapshot.plus=plus;
snapshot.reference=reference;
snapshot.tensor=tensor;
snapshot.ua=ua;
snapshot.ub=ub;
snapshot.reference_error=reference_error;
snapshot.full_order_error=max(abs(answer-full_result.value));
snapshot.full_answer=full_result.value;
snapshot.full_diagnostics=full_diagnostics;

end

% Independent dense reference for the migrated gcce_two_centres example
function [reference_error,snapshot]=gcce_two_centres_reference(answer,diagnostics)

% Magnet and isotopes; zero electron g isolates the model energy scales
sys.magnet=0; sys.isotopes={'E','E','13C','13C'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=repmat({zeros(3)},1,4);
inter.coupling.matrix=cell(4);

% Standard Spinach coupling and zero-field-splitting tensors in Hz
inter.coupling.matrix{1,2}=diag([700 500 1200]);
inter.coupling.matrix{1,3}=diag([1700 900 2000]);
inter.coupling.matrix{2,4}=diag([1300 800 1600]);
inter.coupling.matrix{3,4}=diag([100 100 -200]);

% Compile metadata, without a full-bath basis


spin_system=create(sys,inter);

% Generalised CCE propagates the entire observed density matrix
parameters.observed_spins=[1 2];
parameters.bath_spins=[3 4];
parameters.method='gcce'; parameters.expansion='standard';
parameters.quantity='density'; parameters.order=2;
v=[1;0;0;1]/sqrt(2); parameters.rho0=v*v';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,31); parameters.sequence=struct([]);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Run full order for the small bath, retaining all observed density-matrix off-diagonal elements
result=diagnostics.reduction; result.value=answer;
result.times=parameters.times(:);
result.density=reshape(answer.',sqrt(size(answer,2)),sqrt(size(answer,2)),[]);
result.diagnostics.trace_error=zeros(numel(parameters.times),1);
result.diagnostics.minimum_eigenvalue=zeros(numel(parameters.times),1);

% Check physicality before evaluating spin correlations
negativity=NaN(size(result.times));
for n=1:numel(result.times)

    % Partial transpose of the second observed electron, no eigenvalue clipping
    rho=result.density(:,:,n);
    if all(result.valid(n,:))&&abs(trace(rho)-1)<1e-8&&...
       norm(rho-rho','fro')<1e-8&&min(eig((rho+rho')/2))>-1e-8
        blocks=reshape(rho,2,2,2,2);
        partial=reshape(permute(blocks,[3 2 1 4]),4,4);
        negativity(n)=(sum(abs(eig(partial)))-1)/2;
    end
end
result.negativity=negativity;

% Independently construct the full dense Hamiltonian from the input Hz tensors
dimensions=spin_system.comp.mults; dimension=prod(dimensions);
operators=cell(numel(dimensions),3); h=zeros(dimension);
for j=1:numel(dimensions)

    % Local angular momentum matrices and their tensor-product embeddings
    q=(dimensions(j)-1)/2; m=q:-1:-q;
    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
    for a=1:3
        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
            kron(local{a},eye(prod(dimensions(j+1:end)))));
    end
end
for j=1:numel(dimensions)
    for k=j:numel(dimensions)

        % Spinach ordered pair matrices transpose into Cartesian operators
        tensor=inter.coupling.matrix{j,k};
        if isempty(tensor), continue; end
        for a=1:3
            for b=1:3
                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
            end
        end
    end
end

% Propagate an explicit product density and trace over the whole finite bath
d_obs=size(parameters.rho0,1); db=dimension/d_obs;
rho0=kron(parameters.rho0,eye(db)/db); reference=zeros(size(result.value));
for n=1:numel(result.times)

    % This independent calculation uses expm, not the CCE eigensolver
    u=expm(-1i*h*result.times(n)); rho=u*rho0*u'; reduced=zeros(d_obs);
    for k=1:db
        ids=k:db:dimension; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
result.reference_error=max(abs(result.value-reference),[],'all');

% Plot observed populations without silently repairing the density matrix
populations=zeros(numel(result.times),4);
for n=1:numel(result.times)

    % Read the diagonal in the documented observed tensor-product order
    populations(n,:)=real(diag(result.density(:,:,n))).';
end

% Preserve every frozen physical reference array
reference_error=result.reference_error; snapshot=struct();
snapshot.a=a;
snapshot.b=b;
snapshot.blocks=blocks;
snapshot.d_obs=d_obs;
snapshot.db=db;
snapshot.dimension=dimension;
snapshot.dimensions=dimensions;
snapshot.h=h;
snapshot.ids=ids;
snapshot.local=local;
snapshot.m=m;
snapshot.negativity=negativity;
snapshot.operators=operators;
snapshot.partial=partial;
snapshot.plus=plus;
snapshot.populations=populations;
snapshot.reduced=reduced;
snapshot.reference=reference;
snapshot.rho=rho;
snapshot.rho0=rho0;
snapshot.u=u;
snapshot.v=v;
end


% Independent dense reference for the migrated mecce_dissipative_bath example
function [reference_error,snapshot]=mecce_dissipative_bath_reference(answer,diagnostics)

% Full density ME-gCCE with an explicitly supplied nuclear jump rate
parameters.observed_spins=1; parameters.bath_spins=2;
parameters.method='megcce'; parameters.expansion='standard';
parameters.quantity='density'; parameters.order=1;
v=[1;1]/sqrt(2); parameters.rho0=v*v';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete'); parameters.times=linspace(0,200e-6,31);
parameters.sequence=struct([]); parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);
jump=[0 1;0 0]; rate=3000;
parameters.jumps=struct('support',2,'operator',jump,'rate',rate);

% Propagate the reduced density with the cluster implementation
result=diagnostics.reduction; result.value=answer;
result.times=parameters.times(:);
result.density=reshape(answer.',sqrt(size(answer,2)),sqrt(size(answer,2)),[]);
result.diagnostics.trace_error=zeros(numel(parameters.times),1);
result.diagnostics.minimum_eigenvalue=zeros(numel(parameters.times),1);

% Independently construct the column-vectorised dense master equation
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
h=2*pi*(1100*kron(x,x)+900*kron(y,y)+2200*kron(z,z));
a=kron(eye(2),jump); q=a'*a; unit=eye(4);
liouvillian=-1i*(kron(unit,h)-kron(h.',unit))+...
    rate*(kron(conj(a),a)-(kron(unit,q)+kron(q.',unit))/2);
rho0=kron(parameters.rho0,eye(2)/2); reference=zeros(size(result.value));
for n=1:numel(result.times)

    % Dense exponential followed by an explicit nuclear partial trace
    rho=reshape(expm(liouvillian*result.times(n))*rho0(:),4,4); reduced=zeros(2);
    for k=1:2

        % Tensor order is observed electron, then bath nucleus
        ids=k:2:4; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
result.reference_error=max(abs(result.value-reference),[],'all');

% Display observed signal without fitting a microscopic relaxation time

% Preserve every frozen physical reference array
reference_error=result.reference_error; snapshot=struct();
snapshot.a=a;
snapshot.h=h;
snapshot.ids=ids;
snapshot.jump=jump;
snapshot.liouvillian=liouvillian;
snapshot.rate=rate;
snapshot.reduced=reduced;
snapshot.reference=reference;
snapshot.rho=rho;
snapshot.rho0=rho0;
snapshot.unit=unit;
snapshot.v=v;
snapshot.x=x;
snapshot.y=y;
snapshot.z=z;
end

% Exhaustive density reference for the sampled nuclear bath
function [sampling_error,snapshot]=cce_polarized_bath_reference(answer,diagnostics)

% Standard Spinach magnet, isotope list and illustrative Hz tensors
sys.magnet=0.01; sys.isotopes={'E','1H','1H'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 12000].';
inter.coupling.matrix{2,3}=diag([400 400 -800]);
spin_system=create(sys,inter);

% The same physical populations, averaged exhaustively rather than sampled
parameters=diagnostics.parameters; population=[0.8;0.2];
parameters.cce.bath_state=struct('mode','explicit_product',...
    'rho',{{diag(population),diag(population)}});
[exact,exact_diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');
sampling_error=max(abs(answer-exact));
snapshot=struct('population',population,'exact',exact,...
                'exact_diagnostics',exact_diagnostics,'sampling_error',sampling_error);

end

% Independent analytic transverse-noise reference and original physical arrays
function [reference_error,snapshot]=reference_cce_noise_spectrum(answer)
omega=-spin('13C')*0.01; times=linspace(0,10*2*pi/abs(omega),401);
reference=cos(omega*times(:))/4; x=[0 1;1 0]/2;
reference_error=max(abs(answer-reference));
snapshot=struct('omega',omega,'reference',reference,'x',x,'reference_error',reference_error);
end

% Independent dense reference for the migrated finite-drive example
function [reference_error,snapshot]=gcce_finite_drive_reference(answer,diagnostics)

% Independent two-spin Hamiltonian including the finite rectangular drive
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
v=[1;1]/sqrt(2); rho0=kron(v*v',eye(2)/2);
h=2*pi*(1100*kron(x,x)+900*kron(y,y)+2200*kron(z,z));
drive=4e4*kron(x,eye(2));
times=diagnostics.parameters.times(:); exact=zeros(size(answer));
for n=1:numel(times)

    % Keep the independent expm calculation and its original arithmetic order
    t=times(n);
    duration=min(t,20e-6);
    u=expm(-1i*h*(t-duration))*expm(-1i*(h+4e4*kron(x,eye(2)))*duration);
    fullrho=u*rho0*u'; small=zeros(2);
    for k=1:2
        ids=k:2:4; small=small+fullrho(ids,ids);
    end
    exact(n,:)=small(:).';
end

% Preserve the physical reference arrays retained by the stage-zero capture
reference_error=max(abs(answer-exact),[],'all');
populations=real(answer(:,1:3:end));
snapshot=struct('drive',drive,'duration',duration,'exact',exact,...
    'fullrho',fullrho,'h',h,'ids',ids,'populations',populations,...
    'rho0',rho0,'small',small,'t',t,'u',u,'v',v,'x',x,'y',y,'z',z);

end


