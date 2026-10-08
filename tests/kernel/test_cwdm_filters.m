% Tests substance-local basis settings and direct-sum descriptor invariants.
% Syntax:
%
%                       result=test_cwdm_filters()
%
% Outputs:
%
%     result - regression checks of the per-substance input contract
%
% Two uncoupled three-proton substances provide exact descriptor references.
% Changing any setting in block one must leave block two identical. Separate
% electron-nuclear blocks exercise the three-component IK-DNP depth.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_filters()

% Announce the algebraic target
fprintf('TESTING: CWDM per-substance basis filters (T4)\n');
result=new_test_result('kernel/cwdm_filters','CWDM basis filters',...
                      'Local basis settings cannot change another substance.');

% Build two identical proton substances and one spin-free pool
sys.magnet=0; sys.isotopes=repmat({'1H'},1,6);
inter.chem.parts={[1 2 3],[4 5 6],[]}; inter.chem.concs=[1 0 0];
inter.coupling.scalar=cell(6);
inter.coupling.scalar{1,2}=10; inter.coupling.scalar{4,5}=10;
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};
reference=test_spin_system(sys,inter,bas);
result=test_true(result,'dimensions',isequal(reference.bas.nstates,[64;64;1]),...
                 'each three-proton block has 4^3 states and the pool has one');
result=test_true(result,'offsets',isequal(reference.bas.offsets,[0;64;128;129]),...
                 'offsets are zero-based prefix sums with a terminal dimension');
result=test_true(result,'spin-free descriptor',isequal(size(reference.bas.basis{3}),[1 0]),...
                 'a spin-free substance has one unit row and no spin columns');

% Exercise every local filter and approximation setting
settings={struct('projections',{{0,[],[]}}),...
          struct('longitudinal',{{{'1H'},{},{}}}),...
          struct('zero_quantum',{{{[1 2]},{},{}}}),...
          struct('approximation',{{'IK-0','none','none'}},'inter_level',{{1,[],[]}}),...
          struct('approximation',{{'IK-1','none','none'}},'inter_level',{{2,[],[]}},...
                 'prox_level',{{1,[],[]}},'connectivity',{{'scalar_couplings',[],[]}}),...
          struct('approximation',{{'IK-2','none','none'}},'prox_level',{{1,[],[]}},...
                 'connectivity',{{'full_tensors',[],[]}}),...
          struct('approximation',{{'IK-0','none','none'}},'inter_level',{{1,[],[]}},...
                 'manual',{{logical([1 1 0]),false(0,3),false(0,0)}}),...
          struct('sym_group',{{{'S2'},{},{}}},'sym_spins',{{{[1 2]},{},{}}},...
                 'sym_a1g_only',{{true,true,true}})};
for n=1:numel(settings)
    options=bas; fields=fieldnames(settings{n});
    for k=1:numel(fields)
        options.(fields{k})=settings{n}.(fields{k});
    end
    actual=basis(reference,options);
    result=test_true(result,['independence ' int2str(n)],...
                     isequal(actual.bas.basis{2},reference.bas.basis{2}),...
                     'the untouched substance descriptor is exactly unchanged');
    result=test_true(result,['unit rows ' int2str(n)],...
                     all(cellfun(@(x)nnz(x(1,:))==0,actual.bas.basis)),...
                     'filters retain the unit as the first row of every block');
end

% Compare isotope and global-index longitudinal selection
options=bas; options.longitudinal={{},{'1H'},{}};
by_isotope=basis(reference,options);
options.longitudinal={{},{[4 5 6]},{}};
by_index=basis(reference,options);
result=test_true(result,'longitudinal labels',...
                 isequal(by_isotope.bas.basis,by_index.bas.basis),...
                 'isotope strings and global spin labels select the same local spins');

% Check per-substance IK-DNP correlation depths
sys.isotopes={'E','1H','1H','E','1H','1H'};
inter.coupling.scalar{1,2}=1e6; inter.coupling.scalar{1,3}=1e6;
inter.coupling.scalar{4,5}=1e6; inter.coupling.scalar{4,6}=1e6;
inter.chem.parts={[1 2 3],[4 5 6]}; inter.chem.concs=[1 0];
dnp.formalism='sphten-liouv'; dnp.approximation={'IK-DNP','IK-DNP'};
dnp.inter_level={[1 3 2],[1 3 2]};
reference=test_spin_system(sys,inter,dnp);
dnp.inter_level{1}=[1 2 1]; actual=basis(reference,dnp);
result=test_true(result,'IK-DNP vector',...
                 isequal(actual.bas.basis{2},reference.bas.basis{2})&&...
                 size(actual.bas.basis{1},1)<size(reference.bas.basis{1},1),...
                 'the three-vector changes only its hosting substance');

% Reject scalar broadcasting, wrong cardinalities, and cross-block manual rows
bad={struct('approximation','none'),struct('approximation',{{'none'}}),...
     struct('approximation',{{'none','none'}},'manual',{{false(1,6),false(0,3)}})};
for n=1:numel(bad)
    options=bad{n}; options.formalism='sphten-liouv'; rejected=false;
    try
        basis(reference,options);
    catch err
        rejected=contains(err.message,'one element per chemical substance')||...
                 contains(err.message,'number of columns in bas.manual');
    end
    result=test_true(result,['input rejection ' int2str(n)],rejected,...
                     'each setting has exactly one local entry per substance');
end

% Confirm that membership participates in the descriptor cache identity
options=dnp; options.inter_level={[1 3 2],[1 3 2]};
actual=reference; actual.chem.parts={[4 5 6],[1 2 3]};
actual=basis(actual,options);
result=test_true(result,'membership hash',...
                 ~strcmp(actual.bas.basis_hash,reference.bas.basis_hash),...
                 'equal descriptors with different global membership have different hashes');

end


