% Tests operator-basis construction and expansion helpers. Syntax:
%
%                    result=test_operator_basis_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks single-spin tensor bases, bosonic bases, central
% transitions, Stevens operators, single-transition matrices, expansion
% helpers, and sparse preallocation dimensions.
%
% ilya.kuprov@weizmann.ac.il

function result=test_operator_basis_suite()

% State the operator-basis target of the test
result=new_test_result('kernel/operator_basis_suite',...
                       'Operator-basis construction functions',...
                       'operator-basis helpers must satisfy their defining algebra and reconstruct expanded operators.');

% Irreducible spherical tensors obey [Lz,T(k,m)]=m*T(k,m)
mult=3; L=pauli(mult); T=irr_sph_ten(mult,2); projections=2:-1:-2;
for n=1:numel(T)
    result=test_close(result,['irr_sph_ten commutation m=' int2str(projections(n))],...
                      comm(L.z,T{n}),projections(n)*T{n},1e-13,1e-13,...
                      'irreducible spherical tensors are eigenoperators of commutation with Lz');
end
Tall=irr_sph_ten(mult);
result=test_close(result,'irr_sph_ten all ranks count',numel(Tall),mult^2,0,0,...
                  'a spin with multiplicity m has m^2 single-spin irreducible tensor operators');

% Stevens rank-one zero-projection operator is Lz
result=test_close(result,'stevens rank-one z',stevens(3,1,0),L.z,1e-14,1e-14,...
                  'the rank-one q=0 Stevens operator is the angular momentum z operator');

% Weyl boson operators obey their defining number-operator commutators
W=weyl(4);
result=test_close(result,'weyl c*a=n',W.c*W.a,W.n,1e-14,1e-14,...
                  'creation followed by annihilation gives the population number operator in the truncated basis');
result=test_close(result,'weyl number creation commutator',comm(W.n,W.c),W.c,1e-14,1e-14,...
                  'the creation operator raises the number eigenvalue by one');
result=test_close(result,'weyl number annihilation commutator',comm(W.n,W.a),-W.a,1e-14,1e-14,...
                  'the annihilation operator lowers the number eigenvalue by one');

% Bosonic monomials and orthogonalised monomials must have documented structure
B=boson_mono(3);
result=test_close(result,'boson_mono identity',B{1},speye(3),1e-14,1e-14,...
                  'the first bosonic monomial is C^0 A^0, the unit operator');
result=test_close(result,'boson_mono creation',B{2},weyl(3).c,1e-14,1e-14,...
                  'serpentine state 2 corresponds to C^1 A^0, the creation operator');
Bo=boson_ortho(3);
max_overlap=0;
for n=1:numel(Bo)
    for k=1:(n-1)
        max_overlap=max(max_overlap,abs(hdot(Bo{n},Bo{k})));
    end
end
result=test_close(result,'boson_ortho off-diagonal overlaps',max_overlap,0,1e-12,1e-12,...
                  'Gram-Schmidt bosonic monomials have zero Hilbert-Schmidt mutual overlap');

% Single-transition operators follow serpentine indexing exactly
ST=sin_tran(3);
for n=1:9
    [row,col]=lin2kq(3,n,1);
    result=test_close(result,['sin_tran element ' int2str(n)],ST{n},sparse(row,col,1,3,3),1e-14,1e-14,...
                      'single-transition matrices contain exactly one non-zero element at the serpentine location');
end

% Central-transition operators act only on the central two Zeeman levels
CTz=centrans(4,'z');
CTp=centrans(4,'+');
result=test_close(result,'centrans z',CTz,sparse([2 3],[2 3],[0.5 -0.5],4,4),1e-14,1e-14,...
                  'central-transition z operator has spin-half populations on the middle two levels');
result=test_close(result,'centrans plus',CTp,sparse(2,3,1,4,4),1e-14,1e-14,...
                  'central-transition raising operator connects only the middle transition');

% IST and bosonic expansion helpers must reconstruct their source operators
A=[1 2i;-3 4];
[states,coeffs]=oper2ist(A);
result=test_close(result,'oper2ist reconstruction',ist_reconstruct(2,states,coeffs),A,1e-13,1e-13,...
                  'irreducible spherical tensor expansion must reconstruct the original matrix');
[states,coeffs]=ct2ist(4,'z');
result=test_close(result,'ct2ist reconstruction',ist_reconstruct(4,states,coeffs),centrans(4,'z'),1e-13,1e-13,...
                  'central-transition IST expansion must reconstruct the central-transition operator');
[states,coeffs]=enlev2ist(3,2,'S');
P=zeros(3); P(2,2)=1;
result=test_close(result,'enlev2ist spin projector',ist_reconstruct(3,states,coeffs),P,1e-13,1e-13,...
                  'spin energy-level IST expansion must reconstruct the requested projector');
[states,coeffs]=bos2ist('CA',3);
result=test_close(result,'bos2ist creation-annihilation',ist_reconstruct(3,states,coeffs),weyl(3).c*weyl(3).a,1e-13,1e-13,...
                  'bos2ist must expand the requested creation/annihilation product');
[states,coeffs]=oper2bm(P);
result=test_close(result,'oper2bm reconstruction',bm_reconstruct(3,states,coeffs),P,1e-13,1e-13,...
                  'bosonic monomial expansion must reconstruct the original matrix');
[states,coeffs]=enlev2bm(3,2);
result=test_close(result,'enlev2bm projector',bm_reconstruct(3,states,coeffs),P,1e-13,1e-13,...
                  'bosonic energy-level expansion must reconstruct the requested projector');

% Two-spin IST and sparse preallocation dimensions must match the active basis
sys.magnet=0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
bas.formalism='zeeman-hilb'; bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
T20=twospinist(spin_system,1,2,[2 0],'comm');
T20_ref=sqrt(2/3)*(operator(spin_system,{'Lz','Lz'},{1 2})-...
          (operator(spin_system,{'L+','L-'},{1 2})+operator(spin_system,{'L-','L+'},{1 2}))/4);
result=test_close(result,'twospinist rank two zero',T20,T20_ref,1e-14,1e-14,...
                  'two-spin rank-2 zero-projection IST follows its documented Cartesian product expression');
Apre=mprealloc(spin_system,2);
result=test_close(result,'mprealloc zeeman-hilb size',size(Apre),[4 4],0,0,...
                  'Hilbert-space preallocation dimension is the product of spin multiplicities');
bas.formalism='zeeman-liouv';
spin_system=test_spin_system(sys,inter,bas);
Apre=mprealloc(spin_system,2);
result=test_close(result,'mprealloc zeeman-liouv size',size(Apre),[16 16],0,0,...
                  'Zeeman-Liouville preallocation dimension is the square of the Hilbert dimension');

end

% Reconstruct a matrix from IST coefficients.
function A=ist_reconstruct(mult,states,coeffs)
T=irr_sph_ten(mult);
A=zeros(mult);
for n=1:numel(states)
    A=A+coeffs(n)*T{states(n)+1};
end
end

% Reconstruct a matrix from bosonic monomial coefficients.
function A=bm_reconstruct(nlevels,states,coeffs)
B=boson_mono(nlevels);
A=zeros(nlevels);
for n=1:numel(states)
    A=A+coeffs(n)*B{states(n)+1};
end
end
