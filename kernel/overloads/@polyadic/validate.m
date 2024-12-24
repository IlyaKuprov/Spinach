% Checks the internal structure of a polyadic object and throws an
% error if the object does not meet expectations. Syntax:
%
%                            validate(p)
%
% Parameters:
%
%      p  - a polyadic object
% 
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/validate.m>

function validate(p)

% Check the type
if ~isa(p,'polyadic'), error('this object is not a polyadic.'); end

% Check core array, level 1
if ~iscell(p.cores), error('p.cores must be a cell array.'); end

% Check prefix and suffix arrays, level 1
if ~iscell(p.prefix), error('p.prefix must be a cell array.'); end
if ~iscell(p.suffix), error('p.suffix must be a cell array.'); end

% Check core array, level 2
for n=1:numel(p.cores)
    if ~iscell(p.cores{n})
        error('p.cores structure must be {{numeric,...},...}');
    end
    for k=1:numel(p.cores{n})
        if ~isnumeric(p.cores{n}{k})
            error('p.cores structure must be {{numeric,...},...}');
        end
        if isa(p.cores{n}{k},'polyadic'), validate(p.cores{n}{k}); end
    end
end

% Check prefix and suffix arrays, level 2
for n=1:numel(p.prefix)
    if ~isnumeric(p.prefix{n})
        error('all elements of p.prefix must be numeric.');
    end
    if isa(p.prefix{n},'polyadic'), validate(p.prefix{n}); end
end
for n=1:numel(p.suffix)
    if ~isnumeric(p.suffix{n})
        error('all elements of p.suffix must be numeric.');
    end
    if isa(p.suffix{n},'polyadic'), validate(p.suffix{n}); end
end
    
% Check core dimensions
core_dims=zeros(numel(p.cores),2);
for n=1:numel(p.cores)
    nrows=cellfun(@(x)size(x,1),p.cores{n});
    ncols=cellfun(@(x)size(x,2),p.cores{n});
    core_dims(n,1)=prod(nrows(:));
    core_dims(n,2)=prod(ncols(:));
end
if (~all(core_dims(:,1)==core_dims(1,1)))||...
   (~all(core_dims(:,2)==core_dims(1,2)))
    error('terms in the buffered sum have different dimensions.');
end

% Check prefix and suffix dimensions
if (~isempty(p.prefix))&&(size(p.prefix{end},2)~=core_dims(1,1))
    error('inconsistent dimensions in prefix-core product.');
end
if (~isempty(p.suffix))&&(size(p.suffix{end},1)~=core_dims(1,2))
    error('inconsistent dimensions in core-suffix product.');
end
if numel(p.prefix)>1
    for n=1:(numel(p.prefix)-1)
        if size(p.prefix{n},2)~=size(p.prefix{n+1},1)
            error('inconsistent dimensions in the prefix sequence.');
        end
    end
end
if numel(p.suffix)>1
    for n=1:(numel(p.suffix)-1)
        if size(p.suffix{n},2)~=size(p.suffix{n+1},1)
            error('inconsistent dimensions in the suffix sequence.');
        end
    end
end

% Check the number of terms
if numel(p.cores)>100
    disp('WARNING - a polyadic sum with over 100 terms.');
end
if numel(p.prefix)>100
    disp('WARNING - a prefix sequence with over 100 terms.');
end
if numel(p.suffix)>100
    disp('WARNING - a suffix sequence with over 100 terms.');
end

end

% "Please don't feed the dogs - they do not exist."
%
% A notice on a couple of empty 
% kennels in a street in Moscow.

