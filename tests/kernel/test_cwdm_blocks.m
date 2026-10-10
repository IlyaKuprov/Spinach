% Tests exact direct-sum sparsity of spin generators and detection states.
% Syntax:
%
%                       result=test_cwdm_blocks()
%
% Outputs:
%
%     result - T2 checks for two independent two-spin substances
%
% Independent substances cannot acquire cross-substance coherences under
% Hamiltonian, relaxation, or pulse action. Individual coils must have
% support only in their hosting block. State constructors carry their substance concentrations;
% unweighted detection vectors are constructed with coil_state.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_blocks()

% Announce the exact sparsity target
fprintf('TESTING: CWDM direct-sum block structure (T2)\n');
result=new_test_result('kernel/cwdm_blocks','CWDM block structure',...
                      'Independent substances have exactly disjoint support.');

% Build distinct coupled two-spin substances
sys.magnet=14.1; sys.isotopes={'1H','13C','1H','13C'};
inter.chem.parts={1:2,3:4}; inter.chem.concs=[0.7 0.3];
inter.zeeman.scalar={1,2,3,4}; inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=10; inter.coupling.scalar{3,4}=8;
inter.relaxation={'t1_t2'}; inter.r1_rates={1,2,3,4};
inter.r2_rates={3,4,5,6}; inter.rlx_keep='secular';
inter.equilibrium='zero'; inter.temperature=298;
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
s=assume(test_spin_system(sys,inter,bas),'nmr');

% Check sparse patterns of Hamiltonian, relaxation, and pulse generators
operators={hamiltonian(s),relaxation(s),operator(s,'Lx','1H'),...
           operator(s,'Ly','13C'),operator(s,'Lz','1H','left'),...
           operator(s,'Lz','1H','right')};
labels={'Hamiltonian','relaxation','proton pulse','carbon pulse',...
        'left action','right action'};
for n=1:numel(operators)
    [rows,cols]=find(operators{n});
    row_block=discretize(rows,s.bas.offsets+0.5);
    col_block=discretize(cols,s.bas.offsets+0.5);
    result=test_true(result,labels{n},all(row_block==col_block),...
                     'every stored nonzero lies within one offsets block');
end

% Check both exact and cheap coil constructors on every individual spin
for method={'exact','cheap'}
    for n=1:2
        rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
        outside=setdiff(1:s.bas.offsets(end),rows);
        for k=s.chem.parts{n}
            coil=coil_state(s,'L+',k,method{1});
            result=test_true(result,['coil ' method{1} ' ' int2str(k)],...
                             nnz(coil(outside))==0&&nnz(coil(rows))>0,...
                             'a detection state occupies only its hosting substance');
        end
    end
end

% Identity actions retain their selected block in both sparse formats
for side={'left','right','acomm','comm'}
    scale=1+strcmp(side{1},'acomm');
    if strcmp(side{1},'comm'), scale=0; end
    for format={'csc','xyz'}
        for n=1:2
            rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
            expected=sparse(rows,rows,scale,s.bas.offsets(end),s.bas.offsets(end));
            for label={'E','T0,0'}
                actual=operator(s,label{1},s.chem.parts{n}(1),side{1},format{1});
                if strcmp(format{1},'xyz')
                    actual=sparse(actual(:,1),actual(:,2),actual(:,3),...
                                  s.bas.offsets(end),s.bas.offsets(end));
                end
                result=test_true(result,['identity action ' side{1} ' ' format{1} ' ' label{1} ' ' int2str(n)],...
                                 isequal(actual,expected),'identity acts only in the selected substance');
            end
            actual=operator(s,{'E','E'},num2cell(s.chem.parts{n}),side{1},format{1});
            if strcmp(format{1},'xyz')
                actual=sparse(actual(:,1),actual(:,2),actual(:,3),...
                              s.bas.offsets(end),s.bas.offsets(end));
            end
            result=test_true(result,['identity product action ' side{1} ' ' format{1} ' ' int2str(n)],...
                             isequal(actual,expected),'a local identity product contributes once');
        end
    end
    expected=scale*speye(s.bas.offsets(end));
    result=test_true(result,['identity isotope action ' side{1}],...
                     isequal(operator(s,'E','1H',side{1}),expected),...
                     'one matching spin per substance contributes one local identity');
    result=test_true(result,['identity all action ' side{1}],...
                     isequal(operator(s,'E','all',side{1}),2*expected),...
                     'two matching spins per substance contribute two local identities');
    rows=1:s.bas.offsets(2);
    expected=sparse(rows,rows,2*scale,s.bas.offsets(end),s.bas.offsets(end));
    result=test_true(result,['identity numeric action ' side{1}],...
                     isequal(operator(s,'E',[1 2],side{1}),expected),...
                     'a numeric sum leaves unrelated substances empty');
end

% Identity factors retain their selected block and per-spin multiplicity
for method={'cheap','exact','chem'}
    weights=s.chem.concs;
    for n=1:2
        expected=sparse(s.bas.offsets(n)+1,1,weights(n),s.bas.offsets(end),1);
        for label={'E','T0,0'}
            actual=state(s,label{1},s.chem.parts{n}(1),method{1});
            result=test_true(result,['identity ' method{1} ' ' label{1} ' ' int2str(n)],...
                             isequal(actual,expected),'identity occupies only the selected unit coordinate');
        end
        actual=state(s,{'E','E'},num2cell(s.chem.parts{n}),method{1});
        result=test_true(result,['identity product ' method{1} ' ' int2str(n)],...
                         isequal(actual,expected),'a local identity product contributes once');
    end
    expected=sparse(s.bas.offsets(1:2)+1,1,weights,s.bas.offsets(end),1);
    actual=state(s,'E','1H',method{1});
    result=test_true(result,['identity isotope ' method{1}],isequal(actual,expected),...
                     'an isotope sum contributes once for each matching spin');
    actual=state(s,'E','all',method{1});
    result=test_true(result,['identity all ' method{1}],isequal(actual,2*expected),...
                     'two selected spins in a substance contribute twice its unit');
    actual=state(s,'E',[1 2],method{1});
    expected=sparse(1,1,weights(1),s.bas.offsets(end),1);
    result=test_true(result,['identity numeric ' method{1}],isequal(actual,2*expected),...
                     'a numeric sum leaves unrelated blocks empty');
end

% Level projectors retain the hosting block throughout their tensor expansion
for method={'cheap','exact','chem'}
    for n=1:2
        spins=s.chem.parts{n}; rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
        local_sys=sys; local_sys.isotopes=sys.isotopes(spins);
        local_inter=struct(); local_inter.chem.concs=s.chem.concs(n);
        local_bas=bas; local_bas.approximation={'none'};
        local=test_spin_system(local_sys,local_inter,local_bas);
        expected=sparse(rows,1,state(local,'ZL1',1,method{1}),...
                        s.bas.offsets(end),1);
        actual=state(s,'ZL1',spins(1),method{1});
        result=test_close(result,['level projector ' method{1} ' ' int2str(n)],...
                          actual,expected,1e-14,1e-14,...
                          'identity and non-identity terms equal an independent local state');
        fprintf('CWDM_PROJECTOR %s block=%d error=%.16g\n',method{1},n,norm(actual-expected));
    end
end

% Partner expansion keeps full descriptors while constructing local states
for n=1:2
    spins=s.chem.parts{n};
    [actual,descr]=partner_state(s,{{'L+',spins(1)}},{{{'E','Lz'},spins(2)}});
    for k=1:2
        labels={'E','Lz'}; expected=repmat({'E'},1,4);
        expected(spins)={'L+',labels{k}};
        result=test_true(result,['partner descriptor ' int2str(n) ' ' int2str(k)],...
                         isequal(descr{k},expected),'descriptors retain global spin positions');
        reference=state(s,expected(spins),num2cell(spins));
        result=test_close(result,['partner state ' int2str(n) ' ' int2str(k)],...
                          actual{k},reference,0,0,'each partner combination stays within its substance');
    end
    fprintf('CWDM_PARTNER block=%d states=%d\n',n,numel(actual));
end

% Genuinely cross-substance partner specifications remain invalid
rejected=false;
try
    partner_state(s,{{'L+',1}},{{{'E','Lz'},3}});
catch err
    rejected=strcmp(err.identifier,'Spinach:which_subst:crossSubstance');
end
result=test_true(result,'cross-substance partners',rejected,...
                 'active and partner spins must share a substance');

end
