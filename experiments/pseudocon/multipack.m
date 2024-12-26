% Packs multipole moments from a linear stream into a cell array
% that is arranged by ranks. Syntax:
%
%                     Ilm=multipack(ranks,moments)
%
% Inputs:
%
%      ranks - a vector of spherical ranks present, e.g. [0 1 2]
%
%      moments - a vector of multipole moments for each rank,
%                arranged in a linear stream
%
% Output:
%
%      Ilm   - a cell array of vector corresponding to the multipole
%              moments defined in http://dx.doi.org/10.1039/c6cp05437d
%              
%                for L=0,  one element
%              
%                for L=1,  three elements
%              
%                for L=2,  five elements
%
%              et cetera.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=multipack.m>

function Ilm=multipack(ranks,moments)

% Check consistency
grumble(ranks,moments);

% Set the cell array dimensions
Ilm=cell(size(ranks));

% Unpack the ranks
current_position=0;
for k=1:numel(ranks)
    Ilm{k}=moments((current_position+1):(current_position+2*ranks(k)+1));    
    current_position=current_position+2*ranks(k)+1;
end

end

% Consistency enforcement
function grumble(ranks,moments)
if (~isnumeric(ranks))||(~isreal(ranks))||...
   (~isrow(ranks))||any(mod(ranks,1)~=0)||...
   (numel(unique(ranks))~=numel(ranks))||any(ranks<0)
    error('the elements of ranks must be unique non-negative integers.');
end
if (~isnumeric(moments))
    error('moments must be numeric.');
end
if numel(moments)~=sum(2*ranks+1)
    error('the number of entries in moments does not match the declared ranks.');
end
end

% In a war of ideas it is people who get killed.
%
% Stanislaw Jerzy Lec

