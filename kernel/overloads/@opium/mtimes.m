% Matrix products involving an OPIUM object. Syntax:
%
%                         c=mtimes(a,b)
%
% Parameters:
%
%     a,b   - opia or numerical arrays
%
% Outputs:
%
%     c     - multiplication result
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=opium/mtimes.m>

function c=mtimes(a,b)

% When A is a scalar
if ~isa(a,'opium')&&isnumeric(a)&&isscalar(a)
    
    % Return opium multiplied by A
    c=b; c.coeff=a*c.coeff; return
    
end

% When A is not a scalar
if ~isa(a,'opium')&&isnumeric(a)
    
    % Check dimension
    if size(a,2)~=size(b,1)
        error('matrix dimensions are not consistent.');
    end
    
    % Return A multiplied by opium 
    c=b.coeff*a; return
    
end

% When B is a scalar
if ~isa(b,'opium')&&isnumeric(b)&&isscalar(b)
    
    % Return opium multiplied by A
    c=a; c.coeff=b*c.coeff; return
    
end

% When B is not a scalar
if ~isa(b,'opium')&&isnumeric(b)
    
    % Check dimension
    if size(a,2)~=size(b,1)
        error('matrix dimensions are not consistent.');
    end
    
    % Return A multiplied by opium 
    c=a.coeff*b; return
    
end

% When both are opia
if isa(a,'opium')&&isa(b,'opium')
    
    % Check dimension
    if size(a,2)~=size(b,1)
        error('matrix dimensions are not consistent.');
    end
    
    % Update the coefficient
    c=a; c.coeff=a.coeff*b.coeff; return
    
end

% Complain and bomb out
error('operands must be either numeric or opium objects.');

end

% The Nobel Prize winning paper on Fourier NMR spectroscopy
% by Ernst and Anderson was rejected twice by the Journal of
% Chemical Physics before being accepted by Review of Scien-
% tific Instruments.
%
% https://doi.org/10.1063/1.1719961

