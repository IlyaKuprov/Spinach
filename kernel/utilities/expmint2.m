% Computes the nested matrix exponential double integral:
%
%    Integrate[expm(-i*A*(T-t))*B*
%    Integrate[expm(-i*C*(t-x))*D*expm(-i*E*x),{x,0,t}],{t,0,T}]
%
% This corresponds to the (1,3) block of the exponential of the auxiliary
% matrix (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
%
%               I=expmint2(spin_system,A,B,C,D,E,T)
%
% Parameters:
%
%    A,B,C,D,E - square matrices
%
%    T         - upper limit of the outer integral
%
% Output:
%
%    I         - the integral as above
%    
% aditya.dev@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function I=expmint2(spin_system,A,B,C,D,E,T)

% Check consistency
grumble(A,B,C,D,E,T);

% Zero filler block
Z=sparse(size(A,1),size(A,2));

% Auxiliary matrix
auxmat=[A  -1i*B,     Z; ...
        Z      C  -1i*D; ...
        Z      Z      E];

% Exponentiate the auxiliary matrix
P=propagator(spin_system,auxmat,T);

% Build block extractors
BE1=[speye(size(A))  Z  Z];
BE3=[Z; Z; speye(size(A))];

% Extract
I=BE1*P*BE3;

end

function grumble(A,B,C,D,E,T)
if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||...
   (~isnumeric(D))||(~isnumeric(E))||(~isnumeric(T))
    error('all arguments must be numeric.');
end
if (~ismatrix(A))||(~ismatrix(B))||(~ismatrix(C))||...
   (~ismatrix(D))||(~ismatrix(E))
    error('A, B, C, D, E must be matrices.');
end
if (size(A,1)~=size(A,2))||(size(B,1)~=size(B,2))||...
   (size(C,1)~=size(C,2))||(size(D,1)~=size(D,2))||...
   (size(E,1)~=size(E,2))
    error('all matrices must be square.');
end
if (size(A,1)~=size(B,1))||(size(B,1)~=size(C,1))||...
   (size(C,1)~=size(D,1))||(size(D,1)~=size(E,1))
    error('all matrices must have the same dimension.');
end
if (~isreal(T))||(~isscalar(T))
    error('T must be a real scalar.');
end
end

% Beauty may be in the eye of the beholder, but elegance in 
% equations is best compiled in LaTeX - unless, of course, 
% they're formatted in Word, where elegance yields to exis-
% tential despair.
%
% Anonymous physicist 

