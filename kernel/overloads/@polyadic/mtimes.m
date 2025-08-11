% Performs multiplications involving polyadics. Syntax:
%
%                           C=mtimes(A,B)
%
% Parameters:
%
%    A,B  - a polyadic or a numerical array
%
% Outputs:
%
%      C  - a polyadic or a numerical array
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/mtimes.m>

function C=mtimes(A,B)

% When A is a number
if ~isa(A,'polyadic')&&isnumeric(A)&&isscalar(A)
    
    % Multiply smallest cores in the B buffer
    for n=1:numel(B.cores)
        [~,smallest_core]=min(cellfun(@numel,B.cores{n}));
        B.cores{n}{smallest_core}=A*B.cores{n}{smallest_core};
    end
    C=simplify(B); return

end

% When A is a sparse matrix
if ~isa(A,'polyadic')&&isnumeric(A)&&issparse(A)
    
    % Attach as a prefix to B
    B.prefix=[{A} B.prefix]; C=simplify(B); return
    
end

% When A is a full matrix
if ~isa(A,'polyadic')&&isnumeric(A)
    
    % Issue a recursive call
    C=ctranspose(ctranspose(B)*...
                 ctranspose(A)); return
    
end
    
% When B is a number
if ~isa(B,'polyadic')&&isnumeric(B)&&isscalar(B)

    % Multiply smallest cores in the A buffer
    for n=1:numel(A.cores)
        [~,smallest_core]=min(cellfun(@numel,A.cores{n}));
        A.cores{n}{smallest_core}=A.cores{n}{smallest_core}*B;
    end
    C=simplify(A); return
    
end

% When B is a sparse matrix
if ~isa(B,'polyadic')&&isnumeric(B)&&issparse(B)
    
    % Attach as a suffix to A
    A.suffix=[A.suffix {B}]; C=simplify(A); return
    
end

% When B is a full matrix
if ~isa(B,'polyadic')&&isnumeric(B)
    
    % Multiply by suffixes
    for n=numel(A.suffix):-1:1
        B=A.suffix{n}*B;
    end
    B=full(B);
   
    % Multiply by cores
    C=zeros(size(B));
    for n=1:numel(A.cores)
        C=C+kronm(A.cores{n},B);
    end
    
    % Multiply by prefixes
    for n=numel(A.prefix):-1:1
        C=A.prefix{n}*C;
    end
    C=full(C); return
    
end

% When both are polyadic
if isa(A,'polyadic')&&isa(B,'polyadic')

    % Merge threshold and green light flag
    merge_thresh=1024; can_proceed=true();

    % Does anything get in the way?
    can_proceed=isempty(A.suffix)&&can_proceed;
    can_proceed=isempty(B.prefix)&&can_proceed;

    % Are these single-core polyadics?
    can_proceed=isscalar(A.cores)&&isscalar(B.cores)&&can_proceed;

    % Does the Kronecker product term count match?
    can_proceed=(numel(A.cores{1})==numel(B.cores{1}))&&can_proceed;

    % Deeper check
    if can_proceed

        % Are all dimensions OK?
        for n=1:numel(A.cores{1})
            
            % Inner product compatibility check
            can_proceed=(size(A.cores{1}{n},2)==...
                         size(B.cores{1}{n},1))&&can_proceed;

            % Opia are fine in any case...
            if (~isa(A.cores{1}{n},'opium'))&&...
               (~isa(B.cores{1}{n},'opium'))

                % ...but matrix dimensions must be moderate
                can_proceed=(size(A.cores{1}{n},1)<=merge_thresh)&&can_proceed;
                can_proceed=(size(A.cores{1}{n},2)<=merge_thresh)&&can_proceed;
                can_proceed=(size(B.cores{1}{n},1)<=merge_thresh)&&can_proceed;
                can_proceed=(size(B.cores{1}{n},2)<=merge_thresh)&&can_proceed;

            end

        end

    end
   
    % Do the deed
    if can_proceed
       
        % Term-by-term product of krons
        C=polyadic({cell(1,numel(A.cores{1}))});
        for n=1:numel(A.cores{1})
            C.cores{1}{n}=A.cores{1}{n}*B.cores{1}{n};
        end

        % Inherit prefix and suffix
        if ~isempty(A.prefix), C.prefix=A.prefix; end
        if ~isempty(B.suffix), C.suffix=B.suffix; end

        % Simplify and return
        C=simplify(C); return;

    else
        
        % Place B as a suffix of A
        C=simplify(suffix(A,B)); return

    end
    
end
    
% Complain and bomb out
error('operands must be either numeric or polyadic objects.');
    
end

% "Convention obliges me to self-effacingly declare that any errors or
% omissions are my own. They are not. If you spot any, please alert me
% immediately, so I can work out who to blame."
%
% Preface to a cryptanalysis book

