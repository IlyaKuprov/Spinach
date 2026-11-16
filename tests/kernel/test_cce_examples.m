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
    if ismember(names{n},{'cce_erbium_cawo4','cce_hbn_quadrupoles','cce_methyl_molecular','cce_sic_divacancy','cce_sip_donor'})
        [answer,diagnostics]=feval(names{n});
        switch names{n}
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
    elseif ismember(names{n},{'gcce_sic_clock','gcce_highspin_siv','gcce_votpp'})
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
        end
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


