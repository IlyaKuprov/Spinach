% Random walk on SO(3), isotropic rotational diffusion. Syntax:
%
%                  eulers=rwalk(npts,tau_c,dt)
%
% Parameters:
%
%         npts   - number of points in the trajectory
%
%         tau_c  - isotropic rotational correlation 
%                  time, seconds
%
%         dt     - inter-point spacing, seconds
%
% Outputs:
%
%         eulers - npts x 3 array of Euler angles for
%                  for each trajectory point, radians
%
% Note: the angles are NOT increments relative to the previous
%       points, they are angles relative to the starting point
%       of the trajectory.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rwalk.m>

function eulers=rwalk(npts,tau_c,dt)

% Check consistency
grumble(npts,tau_c,dt);

% Random unit jump sequence
jump_angles=randn(npts,3)/sqrt(3);

% Time step and diffusion coefficient
jump_angles=sqrt(dt/tau_c)*jump_angles;

% Require jump angles to be small
if mean(abs(jump_angles))>(pi/32)
    error('jump angles must be small, reduce your dt.');
end

% Compute DCM trajectory
DCM=zeros(3,3,npts); DCM(:,:,1)=eye(3);
for n=2:npts
    R=[ 0  1  0; -1  0  0;  0  0  0]*jump_angles(n,1)+...
      [ 0  0  1;  0  0  0; -1  0  0]*jump_angles(n,2)+...
      [ 0  0  0;  0  0  1;  0 -1  0]*jump_angles(n,3);
    DCM(:,:,n)=expm(R)*DCM(:,:,n-1);
end

% Preallocate the answer
eulers=zeros(npts,3);

% Compute Euler angle trajectory
parfor n=1:npts
    [alp,bet,gam]=dcm2euler(DCM(:,:,n));
    eulers(n,:)=[alp bet gam];
end

end

% Consistency enforcement
function grumble(npts,tau_c,dt)
if (~isnumeric(npts))||(~isreal(npts))||...
   (~isscalar(npts))||(npts<1)||(mod(npts,1)~=0)
    error('npts must be a positive real integer.');
end
if (~isnumeric(tau_c))||(~isreal(tau_c))||...
   (~isscalar(tau_c))||(tau_c<=0)
    error('tau_c must be a positive real number.');
end
if (~isnumeric(dt))||(~isreal(dt))||...
   (~isscalar(dt))||(dt<=0)
    error('dt must be a positive real number.');
end
end

% You can't help them - they have to
% help themselves.
%
% John "Jigsaw" Kramer

