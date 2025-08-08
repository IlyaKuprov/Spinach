% Simplifies the structure of the polyadic object by reordering buffers,
% dropping inconsequential terms, and flattening nested polyadics where
% possible. Syntax:
%
%                              p=simplify(p)
%
% Parameters:
%
%     p   - a polyadic object
%
% Outputs:
%
%     p   - a polyadic or a numeric object
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/simplify.m>

function p=simplify(p)

% Check consistency
grumble(p);

% Get size information
[nrows,ncols]=size(p);

% Flush trivial polyadics into all-zero sparse matrices
if isempty(p.cores), p=spalloc(nrows,ncols,0); return; end

% Loop until static
changes_made=true();
while changes_made
    
    % Default disposition
    changes_made=false();

    % Simplify prefixes
    for n=1:numel(p.prefix)

        % Recursive call for polyadics
        if isa(p.prefix{n},'polyadic')
            p.prefix{n}=simplify(p.prefix{n});
        end
        
        % A zero prefix flushes into all-zero sparse matrix
        if nnz(p.prefix{n})==0, p=spalloc(nrows,ncols,0); return; end
        
        % Drop unit prefixes
        if iseye(p.prefix{n})
            p.prefix{n}=[]; 
            changes_made=true();
        end
        
        % Absorb opium prefixes
        if isa(p.prefix{n},'opium')
            p=p.prefix{n}.coeff*p;
            p.prefix{n}=[];
            changes_made=true();
        end

        % Absorb scalar prefixes
        if isnumeric(p.prefix{n})&&isscalar(p.prefix{n})
            p=p.prefix{n}*p;
            p.prefix{n}=[];
            changes_made=true();
        end

    end
    p.prefix(cellfun(@isempty,p.prefix))=[];

    % Simplify suffixes
    for n=1:numel(p.suffix)

        % Recursive call for polyadics
        if isa(p.suffix{n},'polyadic')
            p.suffix{n}=simplify(p.suffix{n});
        end
        
        % A zero suffix flushes into all-zero sparse matrix
        if nnz(p.suffix{n})==0, p=spalloc(nrows,ncols,0); return; end
        
        % Drop unit suffixes
        if iseye(p.suffix{n})
            p.suffix{n}=[];
            changes_made=true();
        end
        
        % Absorb opium suffixes
        if isa(p.suffix{n},'opium')
            p=p.suffix{n}.coeff*p;
            p.suffix{n}=[];
            changes_made=true();
        end

        % Absorb scalar suffixes
        if isnumeric(p.suffix{n})&&isscalar(p.suffix{n})
            p=p.suffix{n}*p;
            p.suffix{n}=[];
            changes_made=true();
        end

    end
    p.suffix(cellfun(@isempty,p.suffix))=[];

    % Over sum terms
    for n=1:numel(p.cores)

        % Over Kronecker cores
        for k=1:numel(p.cores{n})
            
            % Recursive call for polyadics
            if isa(p.cores{n}{k},'polyadic')

                % First simplify the polyadic
                p.cores{n}{k}=simplify(p.cores{n}{k});

                % Single-core nested polyadics
                if isempty(p.cores{n}{k}.prefix)&&...
                   isempty(p.cores{n}{k}.suffix)&&...
                   isscalar(p.cores{n}{k}.cores)
                    
                    % Elevate the singleton core
                    p.cores{n}=[p.cores{n}(1:(k-1)) ...
                                p.cores{n}{k}.cores{1} ...
                                p.cores{n}((k+1):end)];
                    changes_made=true();

                end

            end

            % A zero core kills the term
            if nnz(p.cores{n}{k})==0
                p.cores{n}=[]; changes_made=true(); break
            end
            
            % Replace unit cores by opia
            if ~isa(p.cores{n}{k},'opium')&&iseye(p.cores{n}{k})
                p.cores{n}{k}=opium(size(p.cores{n}{k},1),1);
                changes_made=true();
            end

        end
        
    end
    p.cores(cellfun(@isempty,p.cores))=[];
    
    % If no terms left, into all-zero sparse matrix
    if isempty(p.cores), p=spalloc(nrows,ncols,0); return; end

    % Over sum terms
    for n=1:numel(p.cores)

        % Over Kronecker cores
        for k=1:numel(p.cores{n})
            
            % Kron up adjacent opia
            if ((k+1)<=numel(p.cores{n}))&&...
               isa(p.cores{n}{k},'opium')&&...
               isa(p.cores{n}{k+1},'opium')
                p.cores{n}{k}=kron(p.cores{n}{k},p.cores{n}{k+1});
                p.cores{n}(k+1)=[]; changes_made=true(); break
            end
            
        end
        
    end

end

% Matricise single-core polyadics
if isa(p,'polyadic')&&isempty(p.prefix)&&...
   isempty(p.suffix)&&isscalar(p.cores)&&...
   isscalar(p.cores{1})
    p=p.cores{1}{1};
end

end

% Consistency enforcement
function grumble(p)
if ~isa(p,'polyadic')
    error('p must be polyadic.');
end
end

% New scientific findings do not become commonly accepted by
% convincing other scientists, but rather by waiting until
% they have died and a new generation of scientists immedia-
% tely starts with it.
%
% Max Planck

