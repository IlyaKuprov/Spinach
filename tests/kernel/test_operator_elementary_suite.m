% Tests elementary operator generators in kernel/operators. Syntax:
%
%                    result=test_operator_elementary_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks analytic commutation relations, indexing conventions,
% and small explicit matrices for low-level operator constructors.
%
% ilya.kuprov@weizmann.ac.il

function result=test_operator_elementary_suite()

% State the operator target of the test
result=new_test_result('kernel/operator_elementary_suite',...
                       'Elementary operator generators',...
                       'low-level operator constructors must satisfy their defining algebraic identities.');

% Check spin-one angular momentum commutation and ladder definitions
S=pauli(3);
result=test_close(result,'pauli [Sx,Sy]',S.x*S.y-S.y*S.x,1i*S.z,1e-14,1e-14,...
                  'angular momentum components obey [Sx,Sy]=i Sz');
result=test_close(result,'pauli ladder definitions',S.p,S.x+1i*S.y,1e-15,1e-15,...
                  'the raising operator is Sx+i Sy');
result=test_close(result,'pauli unit',S.u,speye(3),1e-15,1e-15,...
                  'the unit operator is the identity matrix in the spin space');

% Check finite-truncation Weyl algebra away from the unavoidable edge state
W=weyl(4);
result=test_close(result,'weyl number identity',W.c*W.a,W.n,1e-14,1e-14,...
                  'creation followed by annihilation gives the population number operator');
result=test_close(result,'weyl creation commutator',W.n*W.c-W.c*W.n,W.c,1e-14,1e-14,...
                  'creation raises boson population by one quantum');
result=test_close(result,'weyl annihilation commutator',W.n*W.a-W.a*W.n,-W.a,1e-14,1e-14,...
                  'annihilation lowers boson population by one quantum');

% Check bosonic monomial serpentine indexing
B=boson_mono(3);
result=test_close(result,'boson_mono identity',B{1},W.u(1:3,1:3),1e-15,1e-15,...
                  'the first bosonic monomial is C^0 A^0');
result=test_close(result,'boson_mono creation',B{2},weyl(3).c,1e-15,1e-15,...
                  'serpentine index two is the creation operator C');
result=test_close(result,'boson_mono annihilation',B{3},weyl(3).a,1e-15,1e-15,...
                  'serpentine index three is the annihilation operator A');

% Check Gram-Schmidt orthogonality without imposing normalisation
B=boson_ortho(3);
gram=zeros(numel(B));
for n=1:numel(B)
    for k=1:numel(B)
        gram(n,k)=trace(full(B{n}'*B{k}));
    end
end
result=test_close(result,'boson_ortho off-diagonal Gram',gram-diag(diag(gram)),zeros(size(gram)),1e-12,1e-12,...
                  'orthogonalised bosonic monomials have zero Hilbert-Schmidt overlap');

% Check single-transition basis indexing from the documented 4x4 map
A=sin_tran(4);
ref=sparse(1,4,1,4,4);
result=test_close(result,'sin_tran index ten',A{10},ref,1e-15,1e-15,...
                  'serpentine index ten is the (1,4) single-transition matrix');
ref=sparse(4,1,1,4,4);
result=test_close(result,'sin_tran index seven',A{7},ref,1e-15,1e-15,...
                  'serpentine index seven is the (4,1) single-transition matrix');

% Check central-transition operators embedded in a spin-3/2 manifold
ct_z=spalloc(4,4,2); ct_z(2,2)=0.5; ct_z(3,3)=-0.5;
ct_p=spalloc(4,4,1); ct_p(2,3)=1;
ct_m=spalloc(4,4,1); ct_m(3,2)=1;
result=test_close(result,'centrans z',centrans(4,'z'),complex(ct_z),1e-15,1e-15,...
                  'central-transition Sz acts only on the middle Zeeman doublet');
result=test_close(result,'centrans plus',centrans(4,'+'),complex(ct_p),1e-15,1e-15,...
                  'central-transition raising maps the lower central state to the upper one');
result=test_close(result,'centrans minus',centrans(4,'-'),complex(ct_m),1e-15,1e-15,...
                  'central-transition lowering maps the upper central state to the lower one');

% Check irreducible spherical tensor projection quantum numbers
T=irr_sph_ten(3,2);
proj=[2 1 0 -1 -2];
L=pauli(3);
for n=1:numel(T)
    result=test_close(result,['IST projection m=' int2str(proj(n))],L.z*T{n}-T{n}*L.z,proj(n)*T{n},1e-12,1e-12,...
                      'irreducible spherical tensors obey [Lz,T(k,m)]=m T(k,m)');
end

% Check the simplest Stevens operator and Hermiticity of observable components
result=test_close(result,'stevens rank zero',stevens(3,0,0),speye(3),1e-15,1e-15,...
                  'the rank-zero Stevens operator is the unit operator');
S_pos=stevens(3,2,+1);
S_neg=stevens(3,2,-1);
result=test_close(result,'stevens positive-q Hermiticity',S_pos,S_pos',1e-14,1e-14,...
                  'positive-q Stevens components are Hermitian observables');
result=test_close(result,'stevens negative-q Hermiticity',S_neg,S_neg',1e-14,1e-14,...
                  'negative-q Stevens components are Hermitian observables');

end


