% Tests angular-momentum and matrix indexing helpers. Syntax:
%
%                    result=test_indexing_roundtrip_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that linear and structured index representations are
% mutually consistent for spherical tensors, Wigner functions, and matrix
% serpentine indexing.
%
% ilya.kuprov@weizmann.ac.il

function result=test_indexing_roundtrip_suite()

% State the indexing target of the test
result=new_test_result('kernel/indexing_roundtrip_suite',...
                       'Indexing conversion functions',...
                       'indexing helpers must be exact inverses on valid integer domains.');

% Linear spin-state indexing is zero-based and ordered by increasing L
I=0:24;
[L,M]=lin2lm(I);
result=test_close(result,'lin2lm/lm2lin inverse',lm2lin(L,M),I,0,0,...
                  'LM state indexing must round-trip exactly');
result=test_close(result,'lin2lm first ranks',L(1:9),[0 1 1 1 2 2 2 2 2],0,0,...
                  'linear LM indexing lists ranks as 0, then the three rank-1 states, then rank 2');
result=test_close(result,'lin2lm first projections',M(1:9),[0 1 0 -1 2 1 0 -1 -2],0,0,...
                  'within each rank, projections are listed in decreasing M order');

% Wigner D-function indexing is one-based and ordered by increasing L, then M, then N
J=1:35;
[Lw,Mw,Nw]=lin2lmn(J);
result=test_close(result,'lin2lmn/lmn2lin inverse',lmn2lin(Lw,Mw,Nw),J,0,0,...
                  'LMN Wigner-function indexing must round-trip exactly');
result=test_close(result,'lin2lmn rank-one block',[Lw(2:10); Mw(2:10); Nw(2:10)],...
                  [ones(1,9); 1 1 1 0 0 0 -1 -1 -1; 1 0 -1 1 0 -1 1 0 -1],0,0,...
                  'rank-one Wigner functions are ordered by decreasing M and, within each M, decreasing N');

% Serpentine matrix indexing has documented triangular scan order
S1=serpentine(4,1);
S0=serpentine(4,0);
S1_ref=[1 3 6 10; 2 5 9 13; 4 8 12 15; 7 11 14 16];
result=test_close(result,'serpentine base one',S1,S1_ref,0,0,...
                  'base-one serpentine indexing follows anti-diagonal triangular ordering');
result=test_close(result,'serpentine base zero',S0,S1_ref-1,0,0,...
                  'base-zero serpentine indexing is the base-one table shifted by one');

% Serpentine k,q coordinates and linear indices must be exact inverses in both bases
N=4;
idx1=1:N^2;
[K1,Q1]=lin2kq(N,idx1,1);
result=test_close(result,'lin2kq/kq2lin base one',kq2lin(N,K1,Q1,1),idx1,0,0,...
                  'base-one serpentine matrix indexing must round-trip exactly');
idx0=0:(N^2-1);
[K0,Q0]=lin2kq(N,idx0,0);
result=test_close(result,'lin2kq/kq2lin base zero',kq2lin(N,K0,Q0,0),idx0,0,0,...
                  'base-zero serpentine matrix indexing must round-trip exactly');

end
