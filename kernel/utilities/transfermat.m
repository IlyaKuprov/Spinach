% Transfer matrix calculation for linear amplifiers. Syntax:
%
%               T=transfermat(amp_inps,amp_outs)
%
% Parameters:
%
%   amp_inps - a matrix with amplifier input vectors as columns
%
%   amp_outs - a matrix with amplifier output vectors as columns
%
% Outputs:
%
%   T        - the transfer matrix, such that amp_outs=T*amp_inps
%              in the least squares sense
% 
% Note: the number of input-output vector pairs should be bigger than
% the number of elements in those vectors.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=transfermat.m>

function T=transfermat(amp_inps,amp_outs)

% Check consistency
grumble(amp_inps,amp_outs);

% Run the SVD pseudoinverse
T=amp_outs/amp_inps;

end

% Consistency enforcement
function grumble(amp_inps,amp_outs)
if (~isnumeric(amp_inps))||(size(amp_inps,2)<size(amp_inps,1))
    error('amp_inps must be a stack of column vectors wider than it is tall.');
end
if (~isnumeric(amp_outs))||(size(amp_outs,2)<size(amp_outs,1))
    error('amp_outs must be a stack of column vectors wider than it is tall.');
end
if size(amp_inps,2)~=size(amp_outs,2)
    error('the number of vectors in amp_inps and amp_outs stacks must be the same.');
end
end

% A good friend will always stab you in the front.
%
% Oscar Wilde

