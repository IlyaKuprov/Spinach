% Distributes an array in the user-specified dimension.
%
% Mathworks, Inc.
% i.kuprov@soton.ac.uk

function A=distrib_dim(A,dim)

% Set the stage
spmd

    % Codistributor with default partitioning
    defpart=codistributor1d.unsetPartition;
    CoD=codistributor1d(dim,defpart,size(A));

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

% The hostages decided amongst themselves that the two to 
% be released would be Hiyech Kanji and Ali-Guil Ghanzafar; 
% the former as she was pregnant and the latter for no other 
% reason than his loud snoring, which kept the other hostages 
% awake at night and irritated the terrorists.
%
% Wikipedia, on the 1980 Iranian Embassy Siege

