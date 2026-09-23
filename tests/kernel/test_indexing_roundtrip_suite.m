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

% Announce the test target
fprintf('TESTING: Indexing conversion functions\n');

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

% Signed integer and sparse inputs return the same numbers in the input class
[L8,M8]=lin2lm(int8(I));
result=test_true(result,'lin2lm int8 class',isa(L8,'int8')&&isa(M8,'int8'),...
                 'integer descriptor blocks must keep their class through lin2lm');
result=test_close(result,'lin2lm int8 values',double([L8; M8]),[L; M],0,0,...
                 'integer inputs must give the same ranks and projections as double inputs');
result=test_true(result,'lm2lin int8 inverse',isequal(lm2lin(L8,M8),int8(I)),...
                 'integer LM indexing must round-trip exactly in the input class');
[L16,M16]=lin2lm(int16(0:255));
result=test_true(result,'lin2lm int16 no saturation',isequal(lm2lin(L16,M16),int16(0:255))&&(min(M16)==-15),...
                 'the integer arithmetic must not saturate the class of the input');
result=test_true(result,'lm2lin int8 top rank',isequal(lm2lin(int8(10),int8(-10)),int8(120)),...
                 'the highest rank whose indices fit the class must be accepted');
[Ls,Ms]=lin2lm(sparse(single([0 1 2 3; 0 0 5 8])));
result=test_true(result,'lin2lm sparse single',issparse(Ls)&&issparse(Ms)&&isa(Ms,'single')&&isequal(full(Ms),[0 1 0 -1; 0 0 1 -2]),...
                 'sparse single descriptors must give sparse single ranks and projections');
try
    lin2lm(uint8(1:3)); unsigned_refused=false;
catch
    unsigned_refused=true;
end
result=test_true(result,'lin2lm refuses unsigned',unsigned_refused,...
                 'projections are signed, so unsigned integer inputs must be refused');
try
    lm2lin(int8(11),int8(0)); overflow_refused=false;
catch
    overflow_refused=true;
end
result=test_true(result,'lm2lin refuses overflow',overflow_refused,...
                 'ranks beyond the last one that the integer class of L holds completely must be refused');
try
    lin2lm(int8(121)); partial_refused=false;
catch
    partial_refused=true;
end
result=test_true(result,'lin2lm refuses partial rank',partial_refused,...
                 'integer indices beyond the last rank that the class holds completely must be refused');
[Lb,Mb]=lin2lm(int64(9007199136250224));
result=test_true(result,'lin2lm root step-back',isequal([Lb Mb],int64([94906264 -94906264])),...
                 'a floating-point root that rounds up to the next integer must be stepped back');

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


