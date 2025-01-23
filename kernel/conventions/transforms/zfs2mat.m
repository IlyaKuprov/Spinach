% Converts D and E zero-field splitting parameters described in
% the abstract of (http://dx.doi.org/10.1063/1.1682294) into a
% spin interaction matrix. Syntax:
%
%                  M=zfs2mat(D,E,alp,bet,gam)
%
% Parameters:
%
%    D,E  - real scalar parameters, Hz
%
%    alp  - alpha Euler angle in radians
%
%    bet  - beta Euler angle in radians
%
%    gam  - gamma Euler angle in radians
%
% Outputs:
%
%     M   - symmetric 3x3 matrix, Hz
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=zfs2mat.m>

function M=zfs2mat(D,E,alp,bet,gam)

% Check consistency
grumble(D,E,alp,bet,gam);

% Compute the matrix in the eigenframe
M=[-D/3+E, 0, 0; 0, -D/3-E, 0; 0, 0, 2*D/3];

% Rotate the molecule
R=euler2dcm(alp,bet,gam); M=R*M*R';

% Tidy up the double precision
M=M-eye(3)*trace(M)/3; M=(M+M')/2;

end

% Consistency enforcement
function grumble(D,E,alp,bet,gam)
if (~isnumeric(D))||(~isreal(D))||(~isscalar(D))||...
   (~isnumeric(E))||(~isreal(E))||(~isscalar(E))||...
   (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('all inputs must be real scalars.');
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

