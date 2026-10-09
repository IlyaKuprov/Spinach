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
% support only in their hosting block. State constructors are unweighted
% at this stage; concentration semantics are tested separately.
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
            coil=state(s,'L+',k,method{1});
            result=test_true(result,['coil ' method{1} ' ' int2str(k)],...
                             nnz(coil(outside))==0&&nnz(coil(rows))>0,...
                             'a detection state occupies only its hosting substance');
        end
    end
end

end


