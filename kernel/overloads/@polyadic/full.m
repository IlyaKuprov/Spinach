% Converts a polyadic representation of a matrix into a full mat-
% rix. Syntax: 
%
%                        answer=full(p)
%
% The function opens up all the Kronecker products and uses full
% arithmetic throughout even if some cores are sparse.
%
% Parameters:
%
%     p - a polyadic object
%
% Outputs:
%
%     answer - a full matrix
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/full.m>

function answer=full(p)

% Process nested polyadics
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        if isa(p.cores{n}{k},'polyadic')
            p.cores{n}{k}=full(p.cores{n}{k});
        end
    end
end
for n=1:numel(p.prefix)
    if isa(p.prefix{n},'polyadic')
        p.prefix{n}=full(p.prefix{n});
    end
end
for n=1:numel(p.suffix)
    if isa(p.suffix{n},'polyadic')
        p.suffix{n}=full(p.suffix{n});
    end
end

% Find matrix dimensions
[nrows,ncols]=size(p);

% Preallocate the answer
answer=zeros(nrows,ncols);

% Loop over the sum
for n=1:numel(p.cores)
    
    % Compute the polyadic
    term=full(p.cores{n}{1});
    for k=2:numel(p.cores{n})
        term=kron(term,p.cores{n}{k});
    end
    answer=answer+term;
    
end

% Multiply by prefixes
for n=numel(p.prefix):-1:1
    answer=p.prefix{n}*answer;
end

% Multiply by suffixes
for n=1:numel(p.suffix)
    answer=answer*p.suffix{n};
end

% Make sure the result is full
answer=full(answer);

end

% Weak men are more likely to be socialists.
%
% http://dx.doi.org/10.1016/j.evolhumbehav.2017.04.001

