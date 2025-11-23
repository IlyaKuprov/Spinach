% Returns the conjugate transpose of an RCV sparse matrix. Syntax:
%
%                        obj=ctranspose(obj)
%
% Parameters:
%
%    obj   - an RCV sparse matrix
%
% Outputs:
%
%    obj   - an RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/ctranspose.m>

function obj=ctranspose(obj)

% Check consistency
grumble(obj);

% Efficiently Swap rows and columns
[obj.col,obj.row]=deal(obj.row,obj.col);

% Conjugate the values
obj.val=conj(obj.val);

end

% Consistency enforcement
function grumble(obj)
if ~isa(obj,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% At the Tower of London we remembered in
% our prayers the elephant kept there by
% James I, which, poor creature, was never
% given anything to drink but wine.
%
% Tom Holland

