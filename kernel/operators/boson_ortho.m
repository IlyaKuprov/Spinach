% Orthogonal bosonic monomials calculated from the bosonic mono-
% mial basis produced by boson_mono(nlevels). Gram-Schmidt or-
% thogonalisation is used without normalisation. Syntax: 
%
%                    B=boson_ortho(nlevels)
% 
% Parameters:
%
%     nlevels - number of bosonic ladder population levels 
%
% Outputs:
%
%        B - a cell array of orthogonal bosonic monomials
%
% ilya.kuprov@weizmann.ac.il
% sarbojoy.das@weizmann.ac.il
% 
% <https://spindynamics.org/wiki/index.php?title=boson_ortho.m>

function B=boson_ortho(nlevels)

% Check consistency
grumble(nlevels);

% Bosonic monomials
B=boson_mono(nlevels);

% Gram-Schmidt
for n=1:numel(B)
    for k=1:(n-1)
        B{n}=B{n}-B{k}*hdot(B{k},B{n})/hdot(B{k},B{k});
    end
end

end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)
    error('nlevels must be a positive integer.');
end
end

% "We're all our own prisons. We are each our own wardens. 
%  We do our own time. Prison Is In Your Mind."
% 
% Charles Manson, to Tench and Ford 
% in Netflix series Mindhunter

