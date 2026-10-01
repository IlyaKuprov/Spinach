% Adds suffix matrices to a polyadic. Anything the polyadic
% multiplies will first be multiplied by the suffix matri-
% ces. Syntax:
%
%                       p=suffix(p,a)
%
% Parameters:
%
%      p   -   polyadic object
%
%      a   -   suffix matrix
%
% Outputs:
%
%      p   -   polyadic object
%
% Note: a suffix can be a polyadic itself.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/suffix.m>

function p=suffix(p,a)

% Check consistency
grumble(p);

% Buffer implicit-core scalar suffix as a dimensioned scaled identity
if ~isa(a,'polyadic')&&isscalar(a)&&any(cellfun(@(term)any(cellfun(@(core)isa(core,'function_handle'),term)),p.cores))
    p.suffix=[p.suffix {opium(size(p,2),a)}]; return
end

% Absorb the suffix
if ~isa(a,'polyadic')&&isscalar(a)
    
    % Multiply the last core
    for n=1:numel(p.cores)
        p.cores{n}{end}=a*p.cores{n}{end};
    end
    
else
    
    % Check the dimensions
    if size(p,2)~=size(a,1)
        error('matrix dimension mismatch.');
    end

    % Update suffix array
    p.suffix=[p.suffix {a}];
    
end

end

% Consistency enforcement
function grumble(p)
if ~isa(p,'polyadic')
    error('p must be polyadic.');
end
end

% All organisations that are not actually right-wing
% will over time become left-wing.
%
% O'Sullivan's First Law


