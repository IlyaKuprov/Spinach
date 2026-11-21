% Test: kernel/clexpn_fixtures | Pinned CCE physical inputs
% Validate the immutable NV physical fixture and its explicit model choices.
% Syntax:
%
%                    result=test_clexpn_fixtures()
%
% Outputs:
%
%     result - regression result with recorded checks and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_clexpn_fixtures()

% Inspect the fixed physical dataset without constructing a spin system
result=new_test_result('kernel/clexpn_fixtures','Pinned CCE physical inputs',...
                       'Source identity, tensor units, original grid, and no hidden projection.');
random_state=rng;
[sys,inter,fixture]=cce_fixture('nv_four_c13','6f7282b6');
result=test_true(result,'Fixed four-carbon system',...
                 isequal(sys.isotopes,{'E3','13C','13C','13C','13C'})&&...
                 sys.magnet==0.05&&isequal(size(fixture.coordinates),[4 3]),...
                 'The physical fixture is not an isotope-abundance realisation.');
result=test_true(result,'Physical inputs only',...
                 isequal(sort(fieldnames(sys)),{'isotopes';'magnet'})&&...
                 ~isfield(inter,'coordinates')&&isequal(rng,random_state),...
                 'No pool settings, random draws, or duplicated coordinate dipoles are introduced.');
result=test_true(result,'Pinned source identity',...
                 strcmp(fixture.source.repository,'MICCoMpy/PyCCE')&&...
                 strcmp(fixture.source.commit,'d74d3cb683bf3caa07231654314e554f440833b4')&&...
                 strcmp(fixture.source.doi,'10.1002/adts.202100254')&&...
                 strcmp(fixture.revision,'6f7282b6'),...
                 'The numerical input source is identified independently of algorithm settings.');
result=test_true(result,'Original total-time grid',...
                 isequal(size(fixture.times),[1 41])&&fixture.times(1)==0&&...
                 fixture.times(end)==0.5e-3&&all(diff(fixture.times)>0),...
                 'All 41 original samples through 0.5 ms remain present.');
result=test_close(result,'Source zero-field splitting',inter.coupling.matrix{1,1},...
                  1000*diag([-960000 -960000 1920000]),0,0,...
                  'The source kHz tensor is converted to Hz exactly once.');
result=test_true(result,'Projection remains explicit',...
                 all(cellfun(@(x)nnz(x(:,1:2))==6,inter.coupling.matrix(1,2:5))),...
                 'Full hyperfine tensors are returned before the example selects observed Sz.');

% Compare every bath dipole against its source constants and displacement
for j=1:4
    for k=j+1:4
        r=(fixture.coordinates(k,:)-fixture.coordinates(j,:))*1e-10;
        u=r/norm(r);
        angular=1e-7*1.054571817e-34*(6.7282853242747445*1e7)^2/norm(r)^3;
        expected=angular*(eye(3)-3*(u.'*u))/(2*pi);
        result=test_close(result,'Pinned bath dipole',inter.coupling.matrix{j+1,k+1},...
                          expected,0,0,'No isotope-table dipole is added to the explicit source tensor.');
    end
end

% Caller edits cannot mutate the next returned physical dataset
original=inter.coupling.matrix{1,2};
inter.coupling.matrix{1,2}(:)=0;
[~,fresh,again]=cce_fixture('nv_four_c13','6f7282b6');
result=test_close(result,'Independent fixture returns',fresh.coupling.matrix{1,2},original,0,0,...
                  'A caller chooses its projection without modifying the source fixture.');
result=test_true(result,'Repeated grid and metadata',isequal(fixture,again),...
                 'Repeated reads retain literal grids and source metadata.');

% Unknown names and revisions must never select an approximate replacement
requests={{'diamond','6f7282b6'},{'nv_four_c13','latest'},...
          {1,'6f7282b6'},{'nv_four_c13',1},{'nv_four_c13',''}};
for n=1:numel(requests)
    rejected=false;
    try
        cce_fixture(requests{n}{1},requests{n}{2});
    catch exception
        rejected=startsWith(exception.identifier,'cce_fixture:');
    end
    result=test_true(result,['Invalid fixture request ' num2str(n)],rejected,...
                     'Only the explicitly pinned physical dataset is accepted.');
end

end


