% Computes the nested matrix exponential integral:
%
% Integrate[expm(-i*A*(T-t))*B*
% Integrate[expm(-i*C*(t-x))*D*expm(-i*E*x),{x,0,t}], {t,0,T}]
%
% This corresponds to the (1,3) block of the exponential of the auxiliary
% matrix, For further info see the paper by Charles van Loan
% (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
%               varargout=expmint2(spin_system,A,B,C,D,E,T)
%
% Parameters:
%    A,C,E   - diagonal blocks
%    B,D     - off-diagonal blocks
%    T       - time
%
% Output:
%    M13     - (1,3) block of the auxiliary matrix exponential
%    or if requested
%   [M11_inv*M13, M11_inv] - where M11_inv is conj transpose of (1, 1)
%   block.
%
% aditya.dev@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function varargout=expmint2(spin_system,A,B,C,D,E,T)

% Check consistency
grumble(A,B,C,D,E,T);

% Zero argument shortcut
if (T==0)||(nnz(B)==0)||(nnz(D)==0)
    M13=spalloc(size(A,1),size(E,2),0);
    if nargout>1
        M11_inv=propagator(spin_system,A,T)';
        varargout{1}=M13;
        varargout{2}=M11_inv;
    else
        varargout{1}=M13;
    end
    return;
end

dimA = size(A,1);
dimC = size(C,1);
dimE = size(E,1);

auxmat = [    A,    -1i*B,    sparse(dimA, dimE); ...
    sparse(dimC, dimA),      C,    -1i*D; ...
    sparse(dimE, dimA),    sparse(dimE, dimC), E];

% Exponentiate the auxiliary matrix
auxmat = propagator(spin_system, auxmat, T);

% construct block extractors
% BE1 extracts 1st block columns (dimA columns)
BE1 = sparse(size(auxmat,1), dimA);
BE1(1:dimA, 1:dimA) = speye(dimA);

% BE3 extracts 3rd block columns (dimE columns)
BE3 = sparse(size(auxmat,1), dimE);
start_col_3 = dimA + dimC + 1;
BE3(start_col_3:(start_col_3+dimE-1), 1:dimE) = speye(dimE);

% Inform the user
report(spin_system, 'processing auxiliary matrix blocks...');

M13=(BE1'*auxmat*BE3);

if nargout>1
    M11_inv=(BE1'*auxmat*BE1)';
    varargout{1}=clean_up(spin_system,M11_inv*M13,spin_system.tols.prop_chop);
    varargout{2}=clean_up(spin_system,M11_inv,spin_system.tols.prop_chop);

    % Reclaim memory
    clear('M11_inv','M13','auxmat','BE1','BE3');
else
    varargout{1}=M13;

    % Reclaim memory
    clear('M13','auxmat','BE1','BE3');
end


end

function grumble(A,B,C,D,E,T)
if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||...
        (~isnumeric(D))||(~isnumeric(E))
    error('All matrix arguments must be numeric.');
end
if (~ismatrix(A))||(~ismatrix(B))||(~ismatrix(C))||...
        (~ismatrix(D))||(~ismatrix(E))
    error('A, B, C, D, E must be matrices..');
end
if size(A,1)~=size(A,2) || size(C,1)~=size(C,2) || size(E,1)~=size(E,2)
    error('A, C, E must be square.');
end
if size(B,1)~=size(A,1) || size(B,2)~=size(C,1)
    error('B must have dimensions matching A and C.');
end
if size(D,1)~=size(C,1) || size(D,2)~=size(E,1)
    error('D must have dimensions matching C and E.');
end
if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))
    error('T must be a real scalar.');
end
end

% Beauty may be in the eye of the beholder, but elegance in equations is
% best compiled in LaTeX -- unless, of course, they're formatted in Word,
% where elegance yields to existential despair.
% --  Anonymous Physicist 

