% Tests cached full-row hashes and spherical-tensor consumers. Syntax:
%
%                       result=test_basis_hash()
%
% Outputs:
%
%     result - hash identity, membership, and truncated consumer checks
%
% Full bases supply independent row-membership and projected-operator
% references for truncated multi-substance bases, including a spin-free pool.
%
% ilya.kuprov@weizmann.ac.il

function result=test_basis_hash()

% Announce the cache contract
fprintf('TESTING: Cached basis row hashes\n');
result=new_test_result('kernel/basis_hash','Cached basis row hashes',...
                      'Full-row hashes preserve local membership and operator coefficients.');

% Build matching full and truncated substance-local bases
sys.magnet=1; sys.isotopes={'1H','2H','1H','2H'};
inter.chem.parts={1:2,3:4,[]}; inter.chem.concs=[1 0 0];
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};
reference=test_spin_system(sys,inter,bas);
bas.approximation={'IK-0','none','none'}; bas.inter_level={1,[],[]};
s=test_spin_system(sys,inter,bas); selection=[];

% Compare stored hashes and membership against numeric reference rows
for n=1:s.bas.nsubst
    descr=s.bas.basis{n}; hashes=s.bas.row_hash{n};
    expected=repmat(' ',size(descr,1),32);
    for k=1:size(descr,1)
        expected(k,:)=md5_hash(full(descr(k,:)));
    end
    result=test_true(result,'stored row identity',isequal(hashes,expected),...
                     'each stored digest hashes exactly the full local descriptor row');
    [present,indices]=hash_ismember(full(descr),hashes);
    result=test_true(result,'dense query membership',all(present)&&...
                     isequal(indices,(1:size(descr,1))'),...
                     'dense queries locate the sparse basis rows in their original order');
    [~,indices]=ismember(descr,reference.bas.basis{n},'rows');
    selection=[selection;reference.bas.offsets(n)+indices]; %#ok<AGROW>
end

% Check missing rows, duplicate references, and empty query/reference arrays
queries=[0 0;1 0;0 1]; hashes=repmat(md5_hash([0 0]),2,1);
[present,indices]=hash_ismember(sparse(queries),hashes);
result=test_true(result,'missing and duplicate rows',...
                 isequal(present,[true;false;false])&&isequal(indices,[1;0;0]),...
                 'duplicates use the first reference row and missing rows return zero');
[present,indices]=hash_ismember(zeros(0,2),hashes);
result=test_true(result,'empty queries',isequal(size(present),[0 1])&&...
                 isequal(size(indices),[0 1]),'empty queries return empty column vectors');
[present,indices]=hash_ismember(queries,repmat(' ',0,32));
result=test_true(result,'empty reference',~any(present)&&~any(indices),...
                 'an empty reference contains no query rows');
[present,indices]=hash_ismember(zeros(2,0),s.bas.row_hash{3});
result=test_true(result,'zero-column rows',all(present)&&all(indices==1),...
                 'spin-free descriptors hash as full one-row zero-column matrices');

% Compare truncated operators with exact full-basis projections
for side={'left','right','comm','acomm'}
    for labels={{'Lx'},{'L+','Lz'}}
        spins=num2cell(1:numel(labels{1}));
        observed=operator(s,labels{1},spins,side{1});
        expected=operator(reference,labels{1},spins,side{1});
        result=test_close(result,'projected operator coefficients',...
                          observed,expected(selection,selection),1e-13,1e-13,...
                          'active-label substitution preserves spectator matching and coefficients');
    end
end

% Compare exact coils with projected full-basis detection states
observed=coil_state(s,{'L+','Lz'},{1,2},'exact');
expected=coil_state(reference,{'L+','Lz'},{1,2},'exact');
result=test_close(result,'truncated coil',observed,expected(selection),1e-13,1e-13,...
                  'a descriptor absent from the truncated basis contributes zero');

% Check reaction lookup with a truncated reactant and complete product
reaction=struct('reactants',1,'products',2,'matching',[1 3;2 4]);
maps=react_gen(s,reaction);
[present,indices]=ismember(s.bas.basis{2},s.bas.basis{1},'rows');
expected=[s.bas.offsets(2)+find(present) indices(present)];
result=test_true(result,'truncated reaction map',isequal(maps{1},expected)&&...
                 any(~present),'missing reactant rows are excluded without changing retained indices');

% Rebuild the compiled context and retain empty hashes through Zeeman conversion
bas.formalism='zeeman-hilb'; bas.approximation={'none','none','none'};
bas=rmfield(bas,'inter_level');
converted=basis(s,bas);
[converted,~,~,~,~]=sim2liouv(converted,struct(),[],[],[]);
result=test_true(result,'Zeeman conversion lifecycle',...
                 all(cellfun(@(x)isequal(size(x),[0 32]),converted.bas.row_hash)),...
                 'rebuilding clears spherical hashes and conversion retains empty Zeeman descriptors');

end


