% Converts D and E zero-field splitting parameters described in
% the abstract of
%
%               http://dx.doi.org/10.1063/1.1682294
%
% into a diagonal spin interaction matrix. Syntax:
%
%                         M=zfs2mat(D,E)
%
% Perameters:
%
%   D,E   - real scalar parameters
%
% Outputs:
%
%     M   - diagonal 3x3 matrix
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=zfs2mat.m>

function M=zfs2mat(D,E)

% Check consistency
grumble(D,E);

% Compute the matrix
M=[-D/3+E, 0, 0; 0, -D/3-E, 0; 0, 0, 2*D/3];

% Tidy up the double precision
M=M-eye(3)*trace(M)/3;

end

% Consistency enforcement
function grumble(D,E)
if (~isnumeric(D))||(~isnumeric(E))
    error('all inputs must be numeric.');
end
if (~isreal(D))||(~isreal(E))
    error('all inputs must be real.');
end
if (numel(D)~=1)||(numel(E)~=1)
    error('all inputs must have one element.');
end
end

% To watch the courageous Afghan freedom fighters battle modern
% arsenals with simple hand-held weapons is an inspiration to 
% those who love freedom. Their courage teaches us a great lesson - 
% that there are things in this world worth defending. To the Afghan
% people, I say on behalf of all Americans that we admire your 
% heroism, your devotion to freedom, and your relentless struggle 
% against your oppressors. 
%
% Ronald Reagan, 21 March 1983

