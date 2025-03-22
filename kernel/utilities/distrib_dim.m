% Distributes an array in the user-specified dimension
% for parallel processing using spmd. Syntax:
%
%                   A=distrib_dim(A,dim)
%
% Parameters:
%
%     A    - a numerical array
%
%     dim  - the distribution dimension
%
% Output:
%
%     A    - a distributed numerical array
%
% Mathworks, Inc.
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=distrib_dim.m>

function A=distrib_dim(A,dim)

% Check consistency
grumble(A,dim);

% Get the size
size_A=size(A);

% Set the stage
spmd

    % Codistributor with default partitioning
    defpart=codistributor1d.unsetPartition;
    CoD=codistributor1d(dim,defpart,size_A);

    % Start local parts
    LocalParts=1;

end

% Retrieve codistribution index ranges
CoD=CoD{1}; partLimits=[0 cumsum(CoD.Partition)];

% Dimension-agnostic index array
idx=repmat({':'},1,ndims(A));

% Pull local parts
for n=1:numel(LocalParts)
    idx{dim}=partLimits(n)+1:partLimits(n+1);
    LocalParts{n}=A(idx{:}); %#ok<AGROW>
end

% Build a distributed array
A=distributed(LocalParts,dim);

end

% Consistency enforcement
function grumble(A,dim)
if ~isnumeric(A), error('A must be numeric.'); end
if (~isreal(dim))||(~isscalar(dim))||...
   (mod(dim,1)~=0)||(dim<1)
    error('dim must be a positive integer.');
end
if dim>ndims(A)
    error('dim exceeds the number of dimensions in A.');
end
end

% The hostages decided amongst themselves that the two to 
% be released would be Hiyech Kanji and Ali-Guil Ghanzafar; 
% the former as she was pregnant and the latter for no other 
% reason than his loud snoring, which kept the other hostages 
% awake at night and irritated the terrorists.
%
% Wikipedia, on the 1980 Iranian Embassy Siege

