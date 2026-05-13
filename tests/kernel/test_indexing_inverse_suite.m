% Tests indexing helper inverses. Syntax:
%
%                    result=test_indexing_inverse_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks serpentine matrix indexing, spin-state L,M indexing,
% and Wigner-function L,M,N indexing over complete low-rank domains.
%
% ilya.kuprov@weizmann.ac.il

function result=test_indexing_inverse_suite()

% State the indexing target of the test
result=new_test_result('kernel/indexing_inverse_suite',...
                       'Indexing inverse helpers',...
                       'indexing helpers must be exact inverses on their finite integer domains.');

% Check the documented base-one serpentine matrix
S_ref=[1 3 6 10;2 5 9 13;4 8 12 15;7 11 14 16];
S=serpentine(4,1);
result=test_close(result,'serpentine base one',S,S_ref,0,0,...
                  'base-one serpentine indexing follows anti-diagonals from lower to upper rows');

% Check the documented base-zero serpentine matrix
S=serpentine(4,0);
result=test_close(result,'serpentine base zero',S,S_ref-1,0,0,...
                  'base-zero serpentine indexing is the base-one table shifted by one');

% Check k,q to linear and back in base-one indexing
[K,Q]=ndgrid(1:4,1:4);
I=kq2lin(4,K,Q,1);
[K_obs,Q_obs]=lin2kq(4,I,1);
result=test_close(result,'kq2lin lin2kq base one K',K_obs,K,0,0,...
                  'every base-one row index must survive a serpentine round-trip');
result=test_close(result,'kq2lin lin2kq base one Q',Q_obs,Q,0,0,...
                  'every base-one column index must survive a serpentine round-trip');
result=test_close(result,'kq2lin base one table',I,S_ref,0,0,...
                  'kq2lin must look up entries in the base-one serpentine table');

% Check k,q to linear and back in base-zero indexing
[K,Q]=ndgrid(0:3,0:3);
I=kq2lin(4,K,Q,0);
[K_obs,Q_obs]=lin2kq(4,I,0);
result=test_close(result,'kq2lin lin2kq base zero K',K_obs,K,0,0,...
                  'every base-zero row index must survive a serpentine round-trip');
result=test_close(result,'kq2lin lin2kq base zero Q',Q_obs,Q,0,0,...
                  'every base-zero column index must survive a serpentine round-trip');
result=test_close(result,'kq2lin base zero table',I,S_ref-1,0,0,...
                  'kq2lin must look up entries in the base-zero serpentine table');

% Build a complete low-rank L,M domain in documented order
L=[]; M=[];
for l=0:5
    for m=l:-1:-l
        L(end+1)=l; %#ok<AGROW>
        M(end+1)=m; %#ok<AGROW>
    end
end

% Check L,M linear indexing and inverse conversion
I=lm2lin(L,M);
[L_obs,M_obs]=lin2lm(I);
result=test_close(result,'lm2lin documented ordering',I,0:(numel(I)-1),0,0,...
                  'linear spin-state indices increase rank by rank and projection by decreasing M');
result=test_close(result,'lm2lin lin2lm L inverse',L_obs,L,0,0,...
                  'every L rank must survive an L,M indexing round-trip');
result=test_close(result,'lm2lin lin2lm M inverse',M_obs,M,0,0,...
                  'every M projection must survive an L,M indexing round-trip');

% Build a complete low-rank L,M,N Wigner-function domain in documented order
L=[]; M=[]; N=[];
for l=0:3
    for m=l:-1:-l
        for n=l:-1:-l
            L(end+1)=l; %#ok<AGROW>
            M(end+1)=m; %#ok<AGROW>
            N(end+1)=n; %#ok<AGROW>
        end
    end
end

% Check L,M,N linear indexing and inverse conversion
I=lmn2lin(L,M,N);
[L_obs,M_obs,N_obs]=lin2lmn(I);
result=test_close(result,'lmn2lin documented ordering',I,1:numel(I),0,0,...
                  'Wigner-function indices start at one and then increase rank by rank');
result=test_close(result,'lmn2lin lin2lmn L inverse',L_obs,L,0,0,...
                  'every Wigner rank must survive an L,M,N indexing round-trip');
result=test_close(result,'lmn2lin lin2lmn M inverse',M_obs,M,0,0,...
                  'every Wigner left projection must survive an L,M,N indexing round-trip');
result=test_close(result,'lmn2lin lin2lmn N inverse',N_obs,N,0,0,...
                  'every Wigner right projection must survive an L,M,N indexing round-trip');

end


