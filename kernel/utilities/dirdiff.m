% Directional derivatives of the matrix exponential. Implements Equation 11
% of Najfeld and Havel (https://doi.org/10.1006/aama.1995.1017) and Equati-
% on 16 of Goodwin and Kuprov (https://doi.org/10.1063/1.4928978). Syntax:
%
%                       D=dirdiff(spin_system,A,B,T,N)
%
% Parameters:
%
%     A - Hamiltonian at the reference point, corresponding
%         to exp(-1i*A*T) propagator
%
%     B - differentiation direction (if a matrix) or direc-
%         tions (if a cell array of matrices)
%
%     T - the time to use in exp(-1i*A*T)
%
%     N - block dimension of the auxiliary matrix, use N=2
%         to get the propagator and its first derivative
%
% Outputs:
%
%     D - a cell array of matrices {D0,D1,D2,...} of Eq 18
%         in Goodwin and Kuprov
%         
% ilya.kuprov@weizmann.ac.uk
% david.goodwin@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=dirdiff.m>

function D=dirdiff(spin_system,A,B,T,N)

% Check consistency
grumble(A,B,T,N);

% Preallocate arrays
auxmat=cell(N,N); D=cell(1,N);

% Build auxiliary matrix
for n=1:N
    for k=1:N
        auxmat{n,k}=sparse(size(A,1),size(A,2));
    end
end
if iscell(B)
    for n=1:(N-1), auxmat{n,n+1}=B{n}; end
else
    for n=1:(N-1), auxmat{n,n+1}=B; end
end
for n=1:N, auxmat{n,n}=A; end

% Tighten up exponentiation tolerance
spin_system.tols.prop_chop=1e-14;

% Exponentiate auxiliary matrix
auxmat=propagator(spin_system,cell2mat(auxmat),T);

% Extract directional derivatives
for n=1:N
    D{n}=factorial(n-1)*auxmat(1:size(A,1),(1:size(A,2))+size(A,2)*(n-1));
end

end

% Consistency enforcement
function grumble(A,B,T,N)
if (~isnumeric(N))||(~isreal(N))||(~isscalar(N))||(N<2)||(mod(N,1)~=0)
    error('N must be a real integer greater than 1.');
end
if iscell(B)
    if numel(B)~=N-1
        error('number of B matrices must equal N-1.');
    else
        for n=1:numel(B)
            if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
                (~isnumeric(B{n}))||(size(B{n},1)~=size(B{n},2))
                error('A and B must be square matrices.');
            end
        end
    end
else
    if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
        (~isnumeric(B))||(size(B,1)~=size(B,2))
        error('A and B must be square matrices.');
    end
end
if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))
    error('T must be a real scalar.');
end
end

% Disciplining yourself to do what you know is right 
% and important, although difficult, is the high road
% to pride, self-esteem, and personal satisfaction.
%
% Margaret Thatcher

