% Roughly averages an interaction tensor with respect to the 
% rotation around a user-specified axis. Syntax:
%
%                       T=axis_tsymm(T,a)
%
% Parameters:
%
%     T - 3x3 real interaction tensor
%
%     a - 3x1 real vector specifying the rotation axis
%
% Output:
%
%     T - 3x3 real interaction tensor
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=axis_tsymm.m>

function A=axis_tsymm(T,a)

% Check consistency
grumble(T,a);

% Preallocate average
A=zeros(3,3);

% Loop over the full rotation
for n=1:360
    R=anax2dcm(a,pi*n/180);
    A=A+R*T*R';
end

% Compute the average
A=A/360;

end

% Consistency enforcement
function grumble(T,a)
if (~isnumeric(a))||(~isreal(a))||...
   (size(a,1)~=3)||(size(a,2)~=1)
	error('d must be a real 3x1 vector.');
end 
if (~isnumeric(T))||(~isreal(T))||...
   (size(T,1)~=3)||(size(T,2)~=3)
	error('T must be a real 3x3 matrix.');
end
end

% A notable American commentator, Charles Krauthammer, once 
% explained Rupert Murdoch's success in founding Fox News,
% a cable channel, by pointing out that he had found a niche
% market - half the country.
%
% The Economist

