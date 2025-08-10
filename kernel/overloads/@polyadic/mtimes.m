% Performs multiplications involving polyadics. Syntax:
%
%                           c=mtimes(a,b)
%
% Parameters:
%
%    a,b  - a polyadic or a numerical array
%
% Outputs:
%
%      c  - a polyadic or a numerical array
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/mtimes.m>

function c=mtimes(a,b)

% When A is a number
if ~isa(a,'polyadic')&&isnumeric(a)&&isscalar(a)
    
    % Multiply smallest cores in the B buffer
    for n=1:numel(b.cores)
        [~,smallest_core]=min(cellfun(@numel,b.cores{n}));
        b.cores{n}{smallest_core}=a*b.cores{n}{smallest_core};
    end
    c=simplify(b); return

end

% When A is a sparse matrix
if ~isa(a,'polyadic')&&isnumeric(a)&&issparse(a)
    
    % Attach as a prefix to B
    b.prefix=[{a} b.prefix]; c=simplify(b); return
    
end

% When A is a full matrix
if ~isa(a,'polyadic')&&isnumeric(a)
    
    % Issue a recursive call
    c=ctranspose(ctranspose(b)*...
                 ctranspose(a)); return
    
end
    
% When B is a number
if ~isa(b,'polyadic')&&isnumeric(b)&&isscalar(b)

    % Multiply smallest cores in the A buffer
    for n=1:numel(a.cores)
        [~,smallest_core]=min(cellfun(@numel,a.cores{n}));
        a.cores{n}{smallest_core}=a.cores{n}{smallest_core}*b;
    end
    c=simplify(a); return
    
end

% When B is a sparse matrix
if ~isa(b,'polyadic')&&isnumeric(b)&&issparse(b)
    
    % Attach as a suffix to A
    a.suffix=[a.suffix {b}]; c=simplify(a); return
    
end

% When B is a full matrix
if ~isa(b,'polyadic')&&isnumeric(b)
    
    % Multiply by suffixes
    for n=numel(a.suffix):-1:1
        b=a.suffix{n}*b;
    end
    b=full(b);
   
    % Multiply by cores
    c=zeros(size(b));
    for n=1:numel(a.cores)
        c=c+kronm(a.cores{n},b);
    end
    
    % Multiply by prefixes
    for n=numel(a.prefix):-1:1
        c=a.prefix{n}*c;
    end
    c=full(c); return
    
end

% When both are polyadic
if isa(a,'polyadic')&&isa(b,'polyadic')
    
    % Add a suffix
    c=simplify(suffix(a,b)); return
    
end
    
% Complain and bomb out
error('operands must be either numeric or polyadic objects.');
    
end

% "Convention obliges me to self-effacingly declare that any errors or
% omissions are my own. They are not. If you spot any, please alert me
% immediately, so I can work out who to blame."
%
% Preface to a cryptanalysis book

