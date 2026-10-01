% Applies a polyadic operator to a numeric vector or matrix. Unlike
% scalar operator scaling, this always returns a numeric block. Syntax:
%
%                         answer=apply(p,x)
%
% Parameters:
%
%    p      - a polyadic operator
%
%    x      - numeric matrix with size(p,2) rows; each column
%             is a separate right-hand side
%
% Outputs:
%
%    answer - numeric matrix with size(p,1) rows and size(x,2) columns
%
% ilya.kuprov@weizmann.ac.il

function answer=apply(p,x)

% Check consistency
grumble(p,x);

% Multiply by suffixes
for n=numel(p.suffix):-1:1
    x=p.suffix{n}*x;
end
x=full(x);

% Preallocate the core product result
cores=core_specs(p);
core_rows=prod(cellfun(@(x)core_size(x,1),cores{1}));
answer=zeros(core_rows,size(x,2),'like',x);

% Multiply by cores
for n=1:numel(p.cores)
    answer=answer+kronm(cores{n},x);
end

% Multiply by prefixes
for n=numel(p.prefix):-1:1
    answer=p.prefix{n}*answer;
end
answer=full(answer);

end

% Consistency enforcement
function grumble(p,x)
if ~isa(p,'polyadic')
    error('p must be polyadic.');
end
if ~isnumeric(x)||~ismatrix(x)||size(x,1)~=size(p,2)
    error('x must be a numeric matrix with size(p,2) rows.');
end
end


