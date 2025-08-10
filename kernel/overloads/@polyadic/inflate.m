% Converts a polyadic representation of a matrix into a sparse mat-
% rix. Syntax: 
%
%                          answer=inflate(p)
%
% The function opens up all the Kronecker products while preserving
% the sparse type if some cores are sparse.
%
% Parameters:
%
%     p - a polyadic object
%
% Outputs:
%
%     answer - a sparse matrix, except if prefixes or suffixes
%              are full (in that case, a full matrix)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/inflate.m>

function answer=inflate(p)

% Process nested polyadics
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        if isa(p.cores{n}{k},'polyadic')
            p.cores{n}{k}=inflate(p.cores{n}{k});
        end
    end
end
for n=1:numel(p.prefix)
    if isa(p.prefix{n},'polyadic')
        p.prefix{n}=inflate(p.prefix{n});
    end
end
for n=1:numel(p.suffix)
    if isa(p.suffix{n},'polyadic')
        p.suffix{n}=inflate(p.suffix{n});
    end
end

% Find matrix dimensions
[nrows,ncols]=size(p);

% Get index arrays going
rows=cell(numel(p.cores),1);
cols=cell(numel(p.cores),1);
vals=cell(numel(p.cores),1);

% Loop over the sum
for n=1:numel(p.cores)
    
    % Compute the polyadic
    term=p.cores{n}{1};
    for k=2:numel(p.cores{n})
        term=kron(term,p.cores{n}{k});
    end
    [rows{n},cols{n},vals{n}]=find(term);
    
end

% Merge index arrays
rows=cell2mat(rows);
cols=cell2mat(cols);
vals=cell2mat(vals);

% Build a sparse matrix
if isempty(vals)
    answer=sparse([],[],[],nrows,ncols);
else
    answer=sparse(rows,cols,vals,nrows,ncols);
end

% Prevent GPU memory leaks
clear('rows','cols','vals');

% Multiply by prefixes
for n=numel(p.prefix):-1:1
    answer=p.prefix{n}*answer;
end

% Multiply by suffixes
for n=1:numel(p.suffix)
    answer=answer*p.suffix{n};
end

end

% Planning to write is not writing. Outlining, researching, 
% talking to people about what you're doing, none of that is
% writing. Writing is writing.
%
% E.L. Doctorow

