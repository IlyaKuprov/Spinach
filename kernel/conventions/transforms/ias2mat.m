% Reconstruction of a 3x3 real interaction matrix C between real 
% vectors u and v from its isotropic-antisymmetric-symmetric de-
% composition:
%
%          a*(u'*v) + d'*cross(u,v) + u'*A*v = u'*C*v
% 
% Syntax:
%
%                      C=ias2mat(a,d,A)
%
% Parameters:
%
%   a - scalar component
%
%   d - antisymmetric coupling vector
%
%   A - symmetric coupling matrix
%
% Outputs:
%
%   C - real 3x3 matrix
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ias2mat.m>

function C=ias2mat(a,d,A)

% Check consistency
grumble(a,d,A);

% Reconstruct the matrix
C=a*eye(3,3)+...
  [ 0     d(3) -d(2);
   -d(3)  0     d(1);
    d(2) -d(1)  0]+A;

end

% Consistency enforcement
function grumble(a,d,A)
if (~isnumeric(a))||(~isreal(a))||(~isscalar(a))
	error('a must be a real scalar.');
end
if (~isnumeric(d))||(~isreal(d))||...
   (size(d,1)~=3)||(size(d,2)~=1)
	error('d must be a real 3x1 vector.');
end 
if (~isnumeric(A))||(~isreal(A))||...
   (size(A,1)~=3)||(size(A,2)~=3)||...
   (norm(A-A',2)/norm(A,2)>1e-6)
	error('A must be a real symmetric 3x3 matrix.');
end 
end

% Мой товарищ, в смертельной агонии
% Не зови понапрасну друзей.
% Дай-ка лучше согрею ладони я
% Над дымящейся кровью твоей.
% 
% Ты не плачь, не стони, ты не маленький,
% Ты не ранен, ты просто убит.
% Я на память сниму с тебя валенки:
% Нам еще наступать предстоит.
% 
% Ион Деген

