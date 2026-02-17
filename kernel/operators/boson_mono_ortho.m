% Orthogonal bosonic monomials calculated from the orginal
% bosonic monomial basis - boson_mono(nlevels). Gram-Schmidt 
% orthogonalization is used to generate B_new from B_old where 
% B_new is an orthogonal basis for bosonic modes. Calculated 
% orthogonal basis has the following property:
%              
%      hdot(B_new{i}, B_new{j})=0, for i not equal to j,
%      hdot(B_new{i}, B_new{j})=N, for i equal to j, N = norm^2 
% where
%                i,j go from 1 to nlevels.
% 
% Syntax: 
%
%                B=boson_mono_ortho(nlevels)
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
% <https://spindynamics.org/wiki/index.php?title=boson_mono_ortho.m>

function B_new = boson_mono_ortho(nlevels)

% Checking consistency of passed argument nlevels
grumble(nlevels);

% Calling the bosonic monomials to be orthogonalized
B_old = boson_mono(nlevels);

% Extracting number of elements
N = numel(B_old); 

% Preallocating cell array for orthogonalized bosonic monomials
B_new = cell(N,1);

% Gram-Schmidt orthogonalization
for n = 1:N
    total_proj = spalloc(nlevels,nlevels,0); 
    for k = 1:(n-1)
        denom = hdot(B_new{k},B_new{k});
        sub_proj = (hdot(B_new{k},B_old{n})/denom)*B_new{k};
        total_proj = total_proj + sub_proj;
    end
        B_new{n} = B_old{n} - total_proj; % Storing orthogonal monomial 
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
% Charles Manson, to Tench and Ford in Netflix series Mindhunter%