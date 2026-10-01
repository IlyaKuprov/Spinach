% Compare implicit Fourier derivatives and polyadic actions with matrices.
% Syntax:
%
%                         result=test_fft_polyadic()
%
% Outputs:
%
%    result - regression checks for FFT actions, adjoints, and propagation
%
% ilya.kuprov@weizmann.ac.il

function result=test_fft_polyadic()

% Describe the independently materialised references
result=new_test_result('kernel/fft_polyadic','FFT polyadic cores',...
                       'Implicit actions must agree with explicit Fourier differentiation and exponentiation.');

% Exercise odd and even grids, including the zero Nyquist convention
for npoints=[1 2 3 8 9]
    [~,explicit]=fourdif(npoints,1);
    core=fourdif_fft(npoints);
    block=reshape(sin(1:(3*npoints)),npoints,3)+...
          1i*reshape(cos(1:(3*npoints)),npoints,3);
    implicit=polyadic({{core,opium(2,1)}});
    reference=kron(explicit,eye(2));
    rhs=[block;block];
    result=test_close(result,['derivative ' num2str(npoints)],...
                      implicit*rhs,reference*rhs,1e-12,1e-12,...
                      'complex multiple right-hand sides match fourdif');
    result=test_close(result,['adjoint ' num2str(npoints)],...
                      implicit'*rhs,reference'*rhs,1e-12,1e-12,...
                      'the adjoint used in norm estimation is correct');
    result=test_close(result,['transpose ' num2str(npoints)],...
                      implicit.'*rhs,reference.'*rhs,1e-12,1e-12,...
                      'non-conjugating transpose has the correct sign');
end

% General complex rectangular handles test dimensions and adjoints
matrix=[1+2i 2-1i 3;4 5i 6-2i];
core=matfree([2 3],@(block)matrix*block,@(block)matrix'*block,false);
implicit=polyadic({{core}});
rhs=[1 2;3 4;5 6];
result=test_close(result,'rectangular action',implicit*rhs,matrix*rhs,...
                  1e-12,1e-12,'forward block dimensions are honoured');
result=test_close(result,'complex adjoint',(2i*implicit)'*[1;2],...
                  (2i*matrix)'*[1;2],1e-12,1e-12,...
                  'scaling conjugates when forming an adjoint');
result=test_close(result,'complex transpose',implicit.'*[1;2],...
                  matrix.'*[1;2],1e-12,1e-12,...
                  'transpose does not conjugate complex entries');
result=test_close(result,'right scaling',(implicit*2)*rhs,2*matrix*rhs,...
                  1e-12,1e-12,'right scalar multiplication remains implicit');
result=test_close(result,'implicit composition',(implicit'*implicit)*rhs,...
                  matrix'*matrix*rhs,1e-12,1e-12,...
                  'composition buffers actions without multiplying cores');
result=test_close(result,'dense left action',[1 2]*implicit,...
                  [1 2]*matrix,1e-12,1e-12,...
                  'dense left multiplication uses the adjoint route');
result=test_true(result,'matrix-free structure',isa(simplify(implicit),'polyadic'),...
                 'a singleton implicit core stays in the polyadic norm-estimation route');

% Explicit materialisation must not open opaque cores
result=check_rejection(result,'full rejection',@()full(implicit),...
                   'matrix-free cores cannot be materialised');
result=check_rejection(result,'inflate rejection',@()inflate(implicit),...
                   'matrix-free cores cannot be materialised');

% Build a small Hermitian rotor generator with a noncommuting spin term
[~,explicit]=fourdif(9,1);
spin_term=[1 0.3;0.3 -1];
reference=kron(diag(cos(2*pi*(0:8)/9)),spin_term)+...
          1i*kron(explicit,eye(2));
implicit=polyadic({{diag(cos(2*pi*(0:8)/9)),spin_term},...
                   {1i*fourdif_fft(9),opium(2,1)}});
spin_system.sys.enable={}; spin_system.sys.disable={};
spin_system.bas.formalism='zeeman-liouv'; spin_system.sys.output='hush';
spin_system.tols.prop_chop=1e-12;
spin_system.tols.dense_matrix=0.5; spin_system.tols.small_matrix=200;
rhs=reshape(sin(1:36),18,2)+1i*reshape(cos(1:36),18,2);
result=test_close(result,'Taylor propagation',...
                  step(spin_system,implicit,rhs,0.1),expm(-1i*reference*0.1)*rhs,...
                  1e-11,1e-11,'the existing step algorithm handles implicit cores');
result=test_close(result,'cleanup preserves action',...
                  clean_up(spin_system,implicit,eps)*rhs,reference*rhs,...
                  1e-12,1e-12,'cleanup does not round opaque actions');

end

% Check that explicit materialisation fails without opening an action
function result=check_rejection(result,label,action,expected)
rejected=false;
try
    action();
catch exception
    rejected=contains(exception.message,expected);
end
result=test_true(result,label,rejected,'opaque cores must not be materialised');
end


