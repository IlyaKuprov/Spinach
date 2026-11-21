% Test: kernel/cce_bath | Physical spin-bath construction
% Deterministic lattice, occupancy, geometry, and frame tests. Syntax:
%
%                    result=test_cce_bath()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_bath()

% Initialise standard harness fields only
result=new_test_result('kernel/cce_bath','CCE bath construction',...
                       'Nonorthogonal cells, occupancy, density and covariance.');

% Compare nonorthogonal replication against direct Cartesian arithmetic
s=struct('mode','lattice','isotopes',{{'1H','13C'}},'magnetic',1,...
         'cell_A',[2 0 0;1 3 0;0.2 0.4 4],'fractional',[0 0 0;0.25 0.5 0.75],...
         'repeats',[2 1 1]);
[sys,~,bath]=spin_bath(struct('magnet',s.magnetic),struct(),rmfield(s,'magnetic'));
expected=[s.fractional*s.cell_A;(s.fractional+[1 0 0])*s.cell_A];
result=test_close(result,'Nonorthogonal row-vector cell convention',bath.coordinates_A,expected,0,0,'Nonorthogonal row-vector cell convention');
result=test_true(result,'Unique physical site IDs',numel(sys.isotopes)==4&&numel(unique(bath.ids))==4,'Unique physical site IDs');

% Check isotope occupancy without mutating the global random stream
s.species={'1H','2H'};
s.occupancy=[0.4 0.6];
s.seed=791;
state=rng;
[~,~,a]=spin_bath(struct('magnet',s.magnetic),struct(),rmfield(s,'magnetic'));
[~,~,b]=spin_bath(struct('magnet',s.magnetic),struct(),rmfield(s,'magnetic'));
result=test_true(result,'Seeded isotope occupancy is reproducible and local',isequal(a,b)&&isequal(state,rng),'Seeded isotope occupancy is reproducible and local');

% Verify active frame rotation of coordinates and Cartesian interactions
q=struct('mode','explicit','isotopes',{{'1H','13C'}},'magnetic',1,...
         'coordinates_A',[1 2 3;4 5 6],'rotation',[0 -1 0;1 0 0;0 0 1],...
         'translation_A',[2 3 4]);
q.pairs_Hz=struct('ids',[1 2],'tensor',[1 2 3;4 5 6;7 8 10]);
[~,inter,bath]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_close(result,'Coordinate covariance',bath.coordinates_A,q.coordinates_A*q.rotation.'+q.translation_A,0,0,'Coordinate covariance');
result=test_close(result,'Tensor covariance',inter.coupling.matrix{1,2},q.rotation*q.pairs_Hz.tensor*q.rotation.',0,0,'Tensor covariance');

% Preserve the preceding geometry fixture while checking tensor conversion
saved_q=q;

% Round-trip local total-frequency tensors through native Zeeman conventions
q=struct('mode','explicit','isotopes',{{'1H','E','15N'}},'magnetic',1,...
         'coordinates_A',[0 0 0;1 0 0;2 0 0]);
q.rotation=[0 -1 0;1 0 0;0 0 1];
q.tensors_Hz=cat(3,[100 2 3;2 200 4;3 4 300],diag([400 500 600]),diag([-70 90 110]));
[sys,inter]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
ss=create(sys,inter);
tensor_error=0;
for k=1:3
    reference=2*pi*q.rotation*q.tensors_Hz(:,:,k)*q.rotation.';
    result=test_close(result,'Local Hz tensor native Zeeman roundtrip',...
                      ss.inter.zeeman.matrix{k},reference,1e-7,1e-10,...
                      'Nuclear ppm and electron g input must reproduce the rotated total angular-frequency tensor.');
    tensor_error=max(tensor_error,norm(ss.inter.zeeman.matrix{k}-reference,'fro'));
end
q.magnetic=0;
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch exception
    rejected=strcmp(exception.identifier,'spin_bath:zero_field');
end
result=test_true(result,'Nonzero local Hz tensor at zero field refused',rejected,...
                 'Native Zeeman inputs cannot represent a field-independent nonzero tensor at zero field.');

% Preserve representable zero tensors at zero magnetic field
q.tensors_Hz=zeros(3,3,3);
[sys,inter]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
ss=create(sys,inter);
result=test_true(result,'Zero local Hz tensors at zero field accepted',...
                 all(cellfun(@(tensor) norm(tensor,'fro')==0,ss.inter.zeeman.matrix)),...
                 'The zero-field refusal is limited to nonzero tensors.');

fprintf('Local Hz roundtrip max rad/s error %.3g.\n',tensor_error);
q=saved_q;

% Select a shell through the documented centre-coordinate field
shell=struct('mode','explicit','isotopes',{{'1H','1H','1H'}},...
             'magnetic',0,'coordinates_A',[0 0 0;1 0 0;2 0 0]);
shell.selection=struct('centres_A',[0 0 0],'inner_A',0.5,'outer_A',1.5);
selected=[];
accepted=false;
try
    [~,~,selected_bath]=spin_bath(struct('magnet',shell.magnetic),struct(),rmfield(shell,'magnetic'));
    selected=selected_bath.ids;
    accepted=true;
catch

    % Convert a validator rejection into a named regression failure
end
result=test_true(result,'Documented shell centres accepted',accepted,...
                 'selection.centres_A must be accepted by the shell validator.');
result=test_true(result,'Documented shell selects physical IDs',isequal(selected,2),...
                 'The shell [0.5,1.5] around the origin must retain only site 2.');

% Test a two-dimensional number density and finite packing exclusion
r=struct('mode','random','isotopes',{{'E'}},'magnetic',0,...
         'bounds_A',[0 0 0;100 200 0],'dimension',2,'density',1e17,...
         'seed',19,'max_attempts',10000,'exclusion_A',1);
[~,~,bath]=spin_bath(struct('magnet',r.magnetic),struct(),rmfield(r,'magnetic'));
surface_density=r.density;
result=test_true(result,'Surface density gives twenty unique sites',numel(bath.ids)==20&&size(unique(bath.coordinates_A,'rows'),1)==20,'Surface density gives twenty unique sites');
result=test_close(result,'SI density and area conversion',bath.realized_density,r.density,1,1e-14,'SI density and area conversion');

% Reject duplicate positions and remove sites before preserving group closure
q=rmfield(q,{'rotation','translation_A','pairs_Hz'});
q.coordinates_A=[0 0 0;1 0 0];
q.groups={[1 2]};
q.remove_ids=2;
[~,~,bath]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_true(result,'Defect removal preserves physical group IDs',isequal(bath.ids,1)&&isequal(bath.groups,{1}),'Defect removal preserves physical group IDs');
q=rmfield(q,'remove_ids');
q.coordinates_A=zeros(2,3);
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch
    rejected=true;
end
result=test_true(result,'Duplicate physical sites are refused',rejected,'Duplicate physical sites are refused');

% Check explicit physical grouping and isotope row normalisation
q=struct('mode','explicit','isotopes',{{'1H','13C'}},'magnetic',1,...
         'coordinates_A',[0 0 0;1 0 0],'groups',{{1,2}});
[sys,inter,b]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_true(result,'Complete groups and isotope row convention',...
                 isequal(b.groups,{1,2})&&isrow(sys.isotopes),...
                 'User groups explicitly partition every retained site.');

% Native input ownership requires an isotope row rather than silent reshaping
malformed=rmfield(q,'magnetic');
malformed.isotopes=malformed.isotopes(:);
refused=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),malformed);
catch
    refused=true;
end
result=test_true(result,'Native isotope column rejected',refused,...
                 'The native bath specification requires explicit isotope row ordering.');

% Check scalar isotropic interactions match explicitly supplied tensors
q.pairs_Hz=struct('ids',[1 2],'tensor',37);
[~,scalar]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
q.pairs_Hz.tensor=37*eye(3);
[~,tensor]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
result=test_close(result,'Scalar isotropic tensor',scalar.coupling.matrix{1,2},...
                  tensor.coupling.matrix{1,2},0,0,'Both inputs are in Hz.');

% Check bath-to-Spinach point-dipole covariance without a Hamiltonian shortcut
q=rmfield(q,'pairs_Hz');
R=[0 -1 0;1 0 0;0 0 1];
q.rotation=R;
q.translation_A=[3 4 5];
[rsys,rinter]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
sys.output='hush';
rsys.output='hush';
sys.parallel={'local',1};
rsys.parallel={'local',1};
plain=create(sys,inter);
rotated=create(rsys,rinter);
T=plain.inter.coupling.matrix{1,2};
rotated_tensor=rotated.inter.coupling.matrix{1,2};
result=test_close(result,'Native dipolar transform covariance',rotated_tensor,R*T*R.',...
                  100*eps*norm(T,'fro'),100*eps,'Rigid translations leave dipoles invariant.');

% Reject invalid original site labels rather than silently dropping them
q=rmfield(q,{'rotation','translation_A'});
q.groups={[1 3]};
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch
    rejected=true;
end
result=test_true(result,'Unknown group site refused',rejected,'Original candidate IDs are checked.');
q.groups={1,2};
q.pairs_Hz=struct('ids',[1 3],'tensor',eye(3));
rejected=false;
try
    spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
catch
    rejected=true;
end
result=test_true(result,'Unknown pair site refused',rejected,'A misspelled interaction cannot disappear.');

% Check the linear density boundary independently of the areal case
r=struct('mode','random','isotopes',{{'1H'}},'magnetic',0,'seed',12,...
         'dimension',1,'bounds_A',[0 2 3;10 2 3],'density',3e9,'max_attempts',20);
[~,~,line]=spin_bath(struct('magnet',r.magnetic),struct(),rmfield(r,'magnetic'));
result=test_true(result,'Linear density and frozen coordinates',numel(line.ids)==3&&...
                 all(line.coordinates_A(:,2)==2)&&all(line.coordinates_A(:,3)==3),...
                 'Three spins in 10 Angstrom at 3e9 per metre.');

% Expose unattainable random packing at the caller ceiling
r.count=2;
r=rmfield(r,'density');
r.exclusion_A=20;
r.max_attempts=2;
rejected=false;
try
    spin_bath(struct('magnet',r.magnetic),struct(),rmfield(r,'magnetic'));
catch e
    rejected=strcmp(e.identifier,'spin_bath:packing');
end
result=test_true(result,'Packing ceiling exposed',rejected,'No hidden relaxation of exclusion radius.');

% Validate joint assignments against count and density-derived random candidates
q=struct('mode','random','isotopes',{{'1H'}},'magnetic',1,'bounds_A',[0 0 0;1 1 1],...
         'count',3,'dimension',3,'seed',17,'max_attempts',100);
q.configurations=struct('isotopes',{{'1H','','13C'},{'','1H','1H'}},'weight',{0.3,0.7});
for density_mode=[false true]
    trial=q;
    if density_mode
        trial=rmfield(trial,'count');
        trial.density=3e30;
    end
    accepted=false;
    try
        [~,~,joint]=spin_bath(struct('magnet',trial.magnetic),struct(),rmfield(trial,'magnetic'));
        accepted=size(joint.candidate_coordinates_A,1)==3&&...
                 isequal(joint.configurations(1).ids,[1 3])&&...
                 isequal(joint.configurations(2).ids,[2 3])&&...
                 isequal([joint.configurations.weight],[0.3 0.7]);
    catch exception
        fprintf('Random joint density=%d: %s\n',density_mode,exception.message);
    end
    result=test_true(result,sprintf('Random joint candidate count density=%d',density_mode),accepted,...
                     'Each joint isotope row covers every candidate before vacancies or exclusions are applied.');

    % Refuse singleton joint rows before any candidate-mask indexing
    trial.configurations=struct('isotopes',{{'1H'}},'weight',1);
    rejected=false;
    try
        spin_bath(struct('magnet',trial.magnetic),struct(),rmfield(trial,'magnetic'));
    catch exception
        rejected=contains(exception.message,'Joint configuration size mismatch');
    end
    result=test_true(result,sprintf('Random joint singleton refused density=%d',density_mode),rejected,...
                     'Malformed joint assignments must fail validation rather than logical-mask indexing.');
end

% Preserve original candidate IDs when random sites are geometrically excluded
q.isotopes={'1H','1H','1H'};
q.remove_ids=2;
accepted=false;
try
    [~,~,joint]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
    accepted=isequal(joint.configurations(1).ids,[1 3])&&...
             isequal(joint.configurations(2).ids,3)&&sum([joint.configurations.weight])==1;
catch exception
    fprintf('Excluded random joint: %s\n',exception.message);
end
result=test_true(result,'Random joint exclusion keeps candidate IDs',accepted,...
                 'Joint rows use the pre-exclusion candidate count even when explicit isotope labels are supplied.');

% Use the same candidate count for site-resolved independent probabilities
q=rmfield(q,{'configurations','remove_ids'});
q.isotopes={'1H'};
q.species={'1H','13C'};
q.occupancy=[1 0;0 1;1 0];
accepted=false;
try
    [sampled,~,~]=spin_bath(struct('magnet',q.magnetic),struct(),rmfield(q,'magnetic'));
    accepted=isequal(sampled.isotopes,{'1H','13C','1H'});
catch exception
    fprintf('Random site probabilities: %s\n',exception.message);
end
result=test_true(result,'Random site-resolved occupancy count',accepted,...
                 'Random per-site independent probabilities have one row per candidate, just like joint configurations.');
fprintf('R5587 count and density each produce three candidates; retained joint IDs [1 3] and [3].\n');

fprintf('CCE bath: lattice error %.3g; density %.16g m^-2; all geometric checks complete.\n',max(abs(expected-a.candidate_coordinates_A),[],'all'),surface_density);

% Preserve absent and explicit Zeeman inputs when appending a tensor-bearing bath
native.isotopes={'E','13C'};
native.magnet=1;
native.output='hush';
native.disable={'hygiene'};
native.parallel={'processes',1};
base=struct('coordinates',{{[0 0 0],[0 0 20]}});
spec=struct('mode','explicit','isotopes',{{'1H'}},'coordinates_A',[0 0 40],...
            'tensors_Hz',diag([100 200 300]));
for variant=1:5
    input=base;
    if variant==2
        input.zeeman=struct();
    end
    if variant==3
        input.zeeman.scalar={2.1,17};
    end
    if variant==4
        input.zeeman.eigs={[2.1 2.2 2.3],[]};
        input.zeeman.euler={[0.2 0.3 0.4],[]};
    end
    if variant==5
        input.zeeman.matrix={diag([2.1 2.2 2.3]),eye(3)*13};
    end
    original=create(native,input);
    [augmented,output]=spin_bath(native,input,spec);
    appended=create(augmented,output);
    for k=1:2
        result=test_close(result,['Appended bath preserves Zeeman ' num2str(variant)],...
                          appended.inter.zeeman.matrix{k},original.inter.zeeman.matrix{k},0,4*eps,...
                          'Introducing a new table preserves existing magnet-frequency or explicit Zeeman tensors.');
    end
end

% Tensor-free appended electrons retain standalone frequencies and dipolar scales
spec=struct('mode','explicit','isotopes',{{'E','13C'}},...
            'coordinates_A',[0 0 40;0 0 60]);
[isolated,isolated_inter]=spin_bath(rmfield(native,'isotopes'),struct(),spec);
standalone=create(isolated,isolated_inter);
for variant=1:7
    input=base;
    switch variant
        case 2
            input.zeeman=struct();
        case 3
            input.zeeman.scalar={2.1,17};
        case 4
            input.zeeman.eigs={[2.1 2.2 2.3],[]};
            input.zeeman.euler={[0.2 0.3 0.4],[]};
        case 5
            input.zeeman.matrix={diag([2.1 2.2 2.3]),eye(3)*13};
        case 6
            input.zeeman.units='Hz';
        case 7
            input.suscept=struct('chi',{{}},'xyz',{{}});
    end
    original=create(native,input);
    [augmented,output]=spin_bath(native,input,spec);
    appended=create(augmented,output);
    for k=1:2
        result=test_close(result,['Tensor-free original Zeeman ' num2str(variant)],...
                          appended.inter.zeeman.matrix{k},original.inter.zeeman.matrix{k},0,4*eps,...
                          'Existing explicit and absent declarations keep their physical meaning.');
        result=test_close(result,['Tensor-free appended Zeeman ' num2str(variant)],...
                          appended.inter.zeeman.matrix{k+2},standalone.inter.zeeman.matrix{k},0,4*eps,...
                          'A tensor-free electron must not lose its magnet-frequency splitting.');
        result=test_close(result,['Tensor-free appended dipolar scale ' num2str(variant)],...
                          appended.inter.zeeman.ddscal{k+2},standalone.inter.zeeman.ddscal{k},0,4*eps,...
                          'Electron magnetic-moment scaling must also retain its standalone value.');
    end
end

% An explicitly supplied zero electron tensor is not replaced by a fallback
spec.isotopes={'E'};
spec.coordinates_A=[0 0 40];
spec.tensors_Hz=zeros(3);
[augmented,output]=spin_bath(native,input,spec);
appended=create(augmented,output);
result=test_true(result,'Explicit appended zero g tensor retained',...
                 nnz(appended.inter.zeeman.matrix{3})==0&&nnz(appended.inter.zeeman.ddscal{3})==0,...
                 'A supplied zero tensor remains distinct from an omitted tensor.');

% Absence of a giant-spin declaration retains the stock no-term default
result=test_true(result,'Absent giant-spin declaration stays absent',...
                 ~isfield(output,'giant')&&all(cellfun(@isempty,appended.inter.giant.coeff)),...
                 'Adding a bath does not introduce a giant-spin model.');

% Append nuclei to an electron with a tilted rank-two giant-spin interaction
native.isotopes={'E3','13C'};
input=base;
input.giant.coeff={{[0 0 0],[0;0;37;0;0]},[]};
input.giant.euler={{[0 0 0],[0.2 0.4 0.7]};{}};
original=create(native,input);
spec=struct('mode','explicit','isotopes',{{'1H'}});
for n=[1 3]
    spec.coordinates_A=[zeros(n,2) (1:n)'*20+20];
    [augmented,output]=spin_bath(native,input,spec);
    result=test_true(result,sprintf('Native giant-spin tables after %d additions',n),...
                     isequal(output.giant.coeff,[input.giant.coeff cell(1,n)])&&...
                     isequal(output.giant.euler,[input.giant.euler;cell(n,1)]),...
                     'Old coefficient and Euler entries retain their values and shapes; new slots are empty.');
    appended=create(augmented,output);
    result=test_true(result,sprintf('Compiled giant-spin terms after %d additions',n),...
                     isequal(appended.inter.giant.coeff,[original.inter.giant.coeff;cell(n,1)]),...
                     'The rotated centre interaction is unchanged and the added nuclei have no giant terms.');
end

end


