% The waveforms on different channels are assumed to be stored in the 
% rows of the input array. The Hessian elements correspond to the ele-
% ments of the waveform array ordered as:
%
%                  [X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn]
%
% where X,Y,Z are different control channels and the index enumerates
% the time discretization points. Gradient dimensions and element or-
% der are the same as the input waveform dimensions and element order.
% Elements of the Hessian are reordered as to correspond to the wavef-
% orm array:
%               [X1 X2 ... Xn Y1 Y2 ... Yn Z1 Z2 ... Zn] 
%
% interchanging the order from controls then time point to time point 
% then controls, or vice versa. Syntax:
%
%                      hess=hess_reorder(hess,K,N)
%
% Parameters:
%
%      hess    -  the old Hessian matrix to be reordered, curre-
%                 ntly ordered K first then N.
%
%      K       -  the first ordered variable of the old Hessian, 
%                 number of control channels in the example above.
%
%      N       -  the second ordered variable of the old Hessian,
%                 number of time points in the example above.
%
% Output:
%
%      hess    -  reordered Hessian with N first then K.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hess_reorder.m>

function hess=hess_reorder(hess,K,N)

% Check consistency
grumble(hess,K,N)

% Do the reordering
hess=reshape(hess,[K N K N]);
hess=permute(hess,[2 1 4 3]);
hess=reshape(hess,[N*K N*K]);

end

% Consistency enforcement
function grumble(hess,dim1,dim2)
if numel(hess)~=(dim1*dim2)^2
    error('Hessian size should be (K*N)^2')
end
end

% If it's true that our species is alone in the universe, 
% then I'd have to say the universe aimed rather low and
% settled for very little.
%
% George Carlin

