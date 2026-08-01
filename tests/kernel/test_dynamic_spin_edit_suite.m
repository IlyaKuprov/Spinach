% Tests deterministic spin-system editing support utilities. Syntax:
%
%                    result=test_dynamic_spin_edit_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks spin removal, dilute-isotope subsystem generation,
% assumption overrides, and merging of Spinach input structures.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_spin_edit_suite()

% Announce the test target
fprintf('TESTING: Spin-system editing utilities\n');

% State the spin editing target of the test
result=new_test_result('kernel/dynamic_spin_edit_suite',...
                       'Spin-system editing utilities',...
                       'local spin-system editing helpers must update dependent metadata without touching external state.');

% Build a small but structurally complete spin-system descriptor
spin_system=local_edit_spin_system();

% Check numeric spin removal updates identities, labels, coordinates, and parts
trimmed=kill_spin(spin_system,2);
result=test_true(result,'kill_spin numeric identities',trimmed.comp.nspins==2&&...
                 isequal(trimmed.comp.isotopes,{'1H','13C'})&&isequal(trimmed.comp.labels,{'h','c2'}),...
                 'removing spin two must delete the second isotope and label, and decrement nspins');
result=test_true(result,'kill_spin numeric dimensions',numel(trimmed.inter.coordinates)==2&&...
                 isequal(size(trimmed.inter.coupling.matrix),[2 2])&&isequal(size(trimmed.inter.proxmatrix),[2 2]),...
                 'coordinate, coupling, and proximity arrays must shrink with the spin count');
result=test_true(result,'kill_spin chemical part update',isequal(trimmed.chem.parts,{[1 2]}),...
                 'chemical part indices must be renumbered after a spin is removed');
result=test_true(result,'kill_spin particle types',isequal(trimmed.comp.types,{'S','S'}),...
                 'particle type codes must shrink with the spin count');
result=test_true(result,'kill_spin isotope hash',strcmp(trimmed.comp.iso_hash,md5_hash({'1H','13C'})),...
                 'the isotope list hash must be recomputed after spin removal');
result=test_true(result,'kill_spin SRSK sources',isequal(trimmed.rlx.srsk_sources,[1 2]),...
                 'scalar relaxation source spins must be renumbered after spin removal');

% Check chemical part renumbering across multiple subsystems
multipart=spin_system; multipart.chem.parts={[1 2],3};
multipart=kill_spin(multipart,3);
result=test_true(result,'kill_spin multi-part update',isequal(multipart.chem.parts,{[1 2],zeros(1,0)}),...
                 'killing the last spin must renumber every chemical subsystem without an error');

% Check destruction of stale basis, connectivity, symmetry, and assumption data
stale=spin_system; stale.bas.formalism='sphten-liouv';
stale.inter.conmatrix=logical(speye(3));
stale.comp.sym_group={'S2'}; stale.comp.sym_spins={[2 3]}; stale.comp.sym_a1g_only=true();
stale.inter.assumptions='nmr';
stale.inter.zeeman.strength={'secular','secular','secular'};
stale.inter.giant.strength={[],[],[]};
stale.inter.coupling.strength=cell(3,3);
stale=kill_spin(stale,2);
result=test_true(result,'kill_spin stale metadata',~isfield(stale,'bas')&&...
                 ~isfield(stale.inter,'conmatrix')&&~isfield(stale.comp,'sym_group')&&...
                 ~isfield(stale.inter,'assumptions')&&~isfield(stale.inter.zeeman,'strength')&&...
                 ~isfield(stale.inter.giant,'strength')&&~isfield(stale.inter.coupling,'strength'),...
                 'basis, connectivity, symmetry, and assumption data must be destroyed on spin removal');

% Check logical spin removal follows the same path
logical_trimmed=kill_spin(spin_system,[false true false]);
result=test_true(result,'kill_spin logical mask',isequal(logical_trimmed.comp.isotopes,trimmed.comp.isotopes)&&...
                 isequal(logical_trimmed.comp.labels,trimmed.comp.labels),...
                 'a logical hit mask must remove the same spin as the equivalent numeric list');

% Check dilute-isotope subsystem generation through kill_spin
subsystems=dilute(spin_system,'13C',1);
carbon_counts=cellfun(@(x)nnz(strcmp(x.comp.isotopes,'13C')),subsystems);
spin_counts=cellfun(@(x)x.comp.nspins,subsystems);
result=test_true(result,'dilute singles',numel(subsystems)==2&&all(carbon_counts==1)&&all(spin_counts==2),...
                 'two labelled 13C spins must produce two two-spin dilute subsystems with one 13C each');

% Check assumption overrides by numeric and isotope specifications
assumed=local_assumed_spin_system();
assumed=dictum(assumed,1,'strong');
assumed=dictum(assumed,[1 2],'weak');
result=test_true(result,'dictum numeric overrides',strcmp(assumed.inter.zeeman.strength{1},'strong')&&...
                 strcmp(assumed.inter.coupling.strength{1,2},'weak'),...
                 'numeric dictum calls must replace Zeeman and coupling strength assumptions in place');
assumed=dictum(assumed,{'13C'},'secular');
result=test_true(result,'dictum isotope override',strcmp(assumed.inter.zeeman.strength{2},'secular'),...
                 'isotope dictum calls must replace the matching Zeeman strength assumptions');

% Check merging of independent Spinach input structure fragments
[sys_parts,inter_parts]=local_merge_parts();
[sys,inter]=merge_inp(sys_parts,inter_parts);
result=test_true(result,'merge_inp sys fields',sys.magnet==14.1&&...
                 isequal(sys.isotopes,{'1H','13C','15N'})&&isequal(sys.labels,{'h','c','n'}),...
                 'merge_inp must concatenate isotope and label row-cell fields while preserving a common magnet');
result=test_true(result,'merge_inp coordinates',numel(inter.coordinates)==3&&...
                 isequal(inter.coordinates{1},[0 0 0])&&isequal(inter.coordinates{3},[0 1 0]),...
                 'merge_inp must concatenate coordinate column-cell fields in subsystem order');
result=test_true(result,'merge_inp common values',isequal(inter.temperature,298)&&...
                 isequal(inter.relaxation,{'redfield'})&&isequal(sys.enable,{'greedy'}),...
                 'merge_inp must keep non-extensive fields that are identical in all subsystems');
result=test_true(result,'merge_inp square blocks',isequal(size(inter.coupling.scalar),[3 3])&&...
                 isequal(inter.coupling.scalar{2,3},5.0)&&isempty(inter.coupling.scalar{1,3})&&...
                 isequal(inter.srfk_mdepth{3,2},7.0)&&isempty(inter.srfk_mdepth{1,2})&&...
                 isequal(inter.weiz_r1d,blkdiag(0.1,[0 0.2; 0.2 0])),...
                 'merge_inp must block-merge square fields with empty cell or zero numeric padding');
result=test_true(result,'merge_inp rate arrays',isequal(inter.r1_rates,{0.1;0.2;0.3})&&...
                 isequal(inter.r2_rates,{1.1;1.2;1.3}),...
                 'merge_inp must merge per-spin rate cell arrays of either orientation into columns');
result=test_true(result,'merge_inp index offsets',isequal(inter.srsk_sources,[1 3])&&...
                 isequal(inter.ignore,{[2 3]})&&isequal(inter.chem.parts,{1,[2 3]})&&...
                 isequal(inter.chem.rates,zeros(2))&&isequal(inter.chem.concs,[1 1]),...
                 'merge_inp must offset spin and subsystem indices by preceding spin counts');
result=test_true(result,'merge_inp suscept centres',isequal(inter.suscept.chi,{0.01*eye(3)})&&...
                 isequal(inter.suscept.xyz,{[5 5 5]}),...
                 'merge_inp must concatenate susceptibility centre lists across subsystems');

% Check column-oriented subsystem lists and partless chemistry
[sys_parts,inter_parts]=local_merge_parts();
inter_parts{2}.chem.parts={1;2};
[~,inter]=merge_inp(sys_parts,inter_parts);
result=test_true(result,'merge_inp column parts',isequal(inter.chem.parts,{1,2,3}),...
                 'column-oriented chemical part lists must merge into offset row lists');
[sys_parts,inter_parts]=local_merge_parts();
inter_parts{1}.chem=struct('rp_theory','haberkorn','rp_electrons',1,'rp_rates',[1e6 2e6]);
inter_parts{2}.chem=struct('rp_theory','haberkorn','rp_electrons',1,'rp_rates',[1e6 2e6]);
inter_parts{1}.tau_c={1e-9}; inter_parts{2}.tau_c={1e-9};
[~,inter]=merge_inp(sys_parts,inter_parts);
result=test_true(result,'merge_inp partless chem',isequal(inter.tau_c,{1e-9})&&...
                 isequal(inter.chem.rp_electrons,[1 2]),...
                 'chemistry without a species split must keep tau_c common and offset electron indices');

% Check that non-extensive differences and malformed inputs are refused
[sys_parts,inter_parts]=local_merge_parts();
inter_parts{2}.temperature=300;
result=test_true(result,'merge_inp value mismatch',local_merge_fails(sys_parts,inter_parts),...
                 'differing non-extensive field values must be refused');
[sys_parts,inter_parts]=local_merge_parts();
inter_parts{2}=rmfield(inter_parts{2},'chem');
result=test_true(result,'merge_inp partial group',local_merge_fails(sys_parts,inter_parts),...
                 'a nested group present in only some subsystems must be refused');
[sys_parts,inter_parts]=local_merge_parts();
inter_parts{1}.bogus=1;
result=test_true(result,'merge_inp unknown field',local_merge_fails(sys_parts,inter_parts),...
                 'unknown inter subfields must be refused');
[sys_parts,inter_parts]=local_merge_parts();
inter_parts{2}.zeeman.bogus=1;
result=test_true(result,'merge_inp unknown group field',local_merge_fails(sys_parts,inter_parts),...
                 'unknown subfields inside nested groups must be refused');
[sys_parts,inter_parts]=local_merge_parts();
sys_parts{2}=rmfield(sys_parts{2},'isotopes');
result=test_true(result,'merge_inp missing isotopes',local_merge_fails(sys_parts,inter_parts),...
                 'sys structures without isotope lists must be refused');

end


% Check that a merge attempt throws an error
function failed=local_merge_fails(sys_parts,inter_parts)
try
    merge_inp(sys_parts,inter_parts); failed=false();
catch
    failed=true();
end
end


function spin_system=local_edit_spin_system()

% Create quiet system output settings
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};

% Define spin identities and basis-independent metadata
spin_system.comp.nspins=3;
spin_system.comp.isotopes={'1H','13C','13C'};
spin_system.comp.types={'S','S','S'};
spin_system.comp.labels={'h','c1','c2'};
spin_system.comp.mults=[2 2 2];
spin_system.comp.iso_hash='unused';

% Define interaction arrays affected by spin removal
spin_system.inter.gammas=[spin('1H') spin('13C') spin('13C')];
spin_system.inter.basefrqs=2*pi*[100e6 25e6 25e6];
spin_system.inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
spin_system.inter.zeeman.ddscal={eye(3),eye(3),eye(3)};
spin_system.inter.giant.coeff={[],[],[]};
spin_system.inter.coupling.matrix=cell(3,3);
spin_system.inter.coordinates={[0 0 0],[1 0 0],[0 1 0]};
spin_system.inter.proxmatrix=zeros(3,3);

% Define relaxation arrays affected by spin removal
spin_system.rlx.r1_rates=[];
spin_system.rlx.r2_rates=[];
spin_system.rlx.lind_r1_rates=[];
spin_system.rlx.lind_r2_rates=[];
spin_system.rlx.srfk_mdepth=[];
spin_system.rlx.weiz_r1d=[];
spin_system.rlx.weiz_r2d=[];
spin_system.rlx.srsk_sources=[1 3];

% Define chemistry arrays affected by spin removal
spin_system.chem.parts={[1 2 3]};
spin_system.chem.flux_rate=[];
spin_system.chem.rp_electrons=[];
spin_system.chem.rp_rates=[];

end


function spin_system=local_assumed_spin_system()

% Start from the edit fixture and add assumption fields
spin_system=local_edit_spin_system();
spin_system.inter.zeeman.strength={'secular','secular','secular'};
spin_system.inter.coupling.strength=cell(3,3);
spin_system.inter.coupling.strength(:)={'secular'};

end


function [sys_parts,inter_parts]=local_merge_parts()

% Build first subsystem input structures
sys_parts{1}.magnet=14.1;
sys_parts{1}.isotopes={'1H'};
sys_parts{1}.labels={'h'};
sys_parts{1}.enable={'greedy'};
inter_parts{1}.coordinates={[0 0 0]};
inter_parts{1}.temperature=298;
inter_parts{1}.relaxation={'redfield'};
inter_parts{1}.zeeman=struct();
inter_parts{1}.coupling.scalar={0.0};
inter_parts{1}.r1_rates={0.1};
inter_parts{1}.r2_rates={1.1};
inter_parts{1}.srfk_mdepth=cell(1,1);
inter_parts{1}.weiz_r1d=0.1;
inter_parts{1}.srsk_sources=1;
inter_parts{1}.ignore={};
inter_parts{1}.suscept.chi={0.01*eye(3)};
inter_parts{1}.suscept.xyz={[5 5 5]};
inter_parts{1}.chem.parts={1};
inter_parts{1}.chem.rates=0;
inter_parts{1}.chem.concs=1;

% Build second subsystem input structures
sys_parts{2}.magnet=14.1;
sys_parts{2}.isotopes={'13C','15N'};
sys_parts{2}.labels={'c','n'};
sys_parts{2}.enable={'greedy'};
inter_parts{2}.coordinates={[1 0 0];[0 1 0]};
inter_parts{2}.temperature=298;
inter_parts{2}.relaxation={'redfield'};
inter_parts{2}.zeeman=struct();
inter_parts{2}.coupling.scalar={0.0 5.0; 0.0 0.0};
inter_parts{2}.r1_rates={0.2;0.3};
inter_parts{2}.r2_rates={1.2;1.3};
inter_parts{2}.srfk_mdepth={[] []; 7.0 []};
inter_parts{2}.weiz_r1d=[0 0.2; 0.2 0];
inter_parts{2}.srsk_sources=2;
inter_parts{2}.ignore={[1 2]};
inter_parts{2}.suscept.chi={};
inter_parts{2}.suscept.xyz={};
inter_parts{2}.chem.parts={[1 2]};
inter_parts{2}.chem.rates=0;
inter_parts{2}.chem.concs=1;

end


