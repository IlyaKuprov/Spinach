% Number of non-zeroes in all kernels of the polyadic. Syntax:
%
%                          answer=nnz(p)
%
% Parameters:
%
%     p - a polyadic object
%
% Outputs:
%
%     answer - an integer number
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/nnz.m>

function answer=nnz(p)

% Start from zero
answer=0;

% Loop over cores
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        answer=answer+nnz(p.cores{n}{k});
    end
end

% Loop over prefix
for n=1:numel(p.prefix)
    answer=answer+nnz(p.prefix{n});
end

% Loop over suffix
for n=1:numel(p.suffix)
    answer=answer+nnz(p.suffix{n});
end

end

% Borderline Probability Disorder: afflicted individuals may 
% dismiss the potential importance of results with P=0.06,
% while unquestioningly accepting the importance of results
% with P=0.05 (see also: significosis).

