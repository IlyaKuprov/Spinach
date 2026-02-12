% Bosonic monomial expansion of specific bosonic energy 
% level projectors. Syntax:
%
%       [states,coeffs]=enlev2bm(nlevels,lvl_num)
%
% Parameters:
%
%     nlevels - number of energy levels in the  mode,
%               a positive integer
%
%     lvl_num - energy level number, counting from the
%               empty mode state upwards
%
% Outputs:
%
%     states  - states, in the Spinach BM basis index-
%               ing, that contribute to the operator in
%               question; use lin2kq to convert to K,Q
%               bosonic monomial indices
%
%     coeffs  - coefficients with which the BMs enter
%               the linear combination
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=enlev2bm.m>

function [states,coeffs]=enlev2bm(nlevels,lvl_num)

% Check consistency
grumble(nlevels,lvl_num);

% Energy level projector
P=zeros(nlevels,nlevels); 
P(lvl_num,lvl_num)=1;

% Bosonic monomial expansion
[states,coeffs]=oper2bm(P);

end

% Consistency enforcement
function grumble(nlevels,lvl_num)
if (~isnumeric(nlevels))||(~isscalar(nlevels))||...
   (~isreal(nlevels))||(nlevels<1)
    error('nlevels must be a positive integer.');
end
if (~isnumeric(lvl_num))||(~isscalar(lvl_num))||...
   (~isreal(lvl_num))||(lvl_num<1)||(lvl_num>nlevels)
    error('the specified bosonic energy level does not exist.');
end
end

% Истребление зануд - долг каждого порядочного 
% человека. Если зануда не разъярён - это позор 
% для окружающих.
%
% Лев Ландау

