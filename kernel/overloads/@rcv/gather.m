% Gathers an RCV sparse matrix from GPU. Syntax:
%
%                    A=gather(A)
%
% Parameters:
%
%    A   - an RCV sparse matrix
%
% Outputs:
%
%    A   - the same matrix with data stored on the CPU
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/gather.m>

function A=gather(A)

% Check consistency
grumble(A);

% Gather to CPU
if A.isGPU
    A.row=gather(A.row);
    A.col=gather(A.col);
    A.val=gather(A.val);
    A.isGPU=false;
end

end

% Consistency enforcement
function grumble(A)
if ~isa(A,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% Aerie, I've noticed the unfortunate fact that you live
% by one of the great lessons of history that nothing is
% often a good thing to do and a clever thing to say.
%
% Edwin Odesseiron, in Baldur's Gate 2

