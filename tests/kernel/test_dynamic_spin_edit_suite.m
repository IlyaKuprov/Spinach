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

end


function spin_system=local_edit_spin_system()

% Create quiet system output settings
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};

% Define spin identities and basis-independent metadata
spin_system.comp.nspins=3;
spin_system.comp.isotopes={'1H','13C','13C'};
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
inter_parts{1}.coordinates={[0 0 0]};
inter_parts{1}.zeeman=struct();
inter_parts{1}.coupling=struct();

% Build second subsystem input structures
sys_parts{2}.magnet=14.1;
sys_parts{2}.isotopes={'13C','15N'};
sys_parts{2}.labels={'c','n'};
inter_parts{2}.coordinates={[1 0 0];[0 1 0]};
inter_parts{2}.zeeman=struct();
inter_parts{2}.coupling=struct();

end


