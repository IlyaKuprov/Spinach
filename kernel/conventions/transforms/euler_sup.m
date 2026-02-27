% Superposition of ZYZ active Euler rotations. Syntax:
%
%          rot_cmp=euler_sup(rot_one,rot_two)
%
% Parameters:
%
%    rot_one   - first Euler angle set [alpha beta gamma],
%                radians, ZYZ active convention
%
%    rot_two   - second Euler angle set [alpha beta gamma],
%                radians, ZYZ active convention
%
% Outputs:
%
%    rot_cmp   - row vector [alpha beta gamma] of the
%                composite rotation, radians
%
% Note: rotations are applied in the supplied order
%
%                 v_rot=R_two*R_one*v
%
%       therefore the composite matrix is
%
%                 R_comp=R_two*R_one
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=euler_sup.m>

function rot_cmp=euler_sup(rot_one,rot_two)

% Check consistency
grumble(rot_one,rot_two);

% Cast to row vectors
rot_one=rot_one(:)';
rot_two=rot_two(:)';

% Wrap both angle sets
rot_one=wrapToPi(rot_one);
rot_two=wrapToPi(rot_two);

% Build DCMs for both
dcm_one=euler2dcm(rot_one);
dcm_two=euler2dcm(rot_two);

% Compose the rotations
dcm_comp=dcm_two*dcm_one;

% Handle identity rotation special case
if all(abs(diag(dcm_comp)-[1;1;1])<1e-12)

    % Identity rotation
    rot_cmp=[pi/2 0 -pi/2];

% Handle the singular beta=pi branch 
elseif (abs(dcm_comp(3,3)+1)<1e-12)&&...
       (abs(rot_one(1)-rot_two(1))<1e-12)&&...
       (abs(rot_one(3)-rot_two(3))<1e-12)

    % Clamp the cosine argument 
    dcm_11=max(min(dcm_comp(1,1),1),-1);

    % Recover the third angle 
    gam=(acos(dcm_11)-pi)/2;

    % Return special case
    rot_cmp=[-gam pi gam];

else

    % Process the rest normally
    rot_cmp=dcm2euler(dcm_comp);

end

% Wrap output into (-pi,pi]
rot_cmp=wrapToPi(rot_cmp);

end

% Consistency enforcement
function grumble(rot_one,rot_two)
if (~isnumeric(rot_one))||(~isreal(rot_one))||(~isvector(rot_one))||...
   (numel(rot_one)~=3)||any(~isfinite(rot_one(:)))
    error('rot_one must be a real finite vector with three elements.');
end
if (~isnumeric(rot_two))||(~isreal(rot_two))||(~isvector(rot_two))||...
   (numel(rot_two)~=3)||any(~isfinite(rot_two(:)))
    error('rot_two must be a real finite vector with three elements.');
end
end

% "My mathematics work is proceeding beyond my wildest
%  hopes, and I am even a bit worried. If it's only in
%  prison that I work so well, will I have to arrange
%  to spend two or three months locked up every year?"
%
% André Weil, in a 1940 letter to his wife,
% after he was first imprisoned in Finland
% (accused of spying) and then, after being
% returned to France, put into prison again
% (convicted of being a deserter).
%
% "We're not all lucky enough to sit and work
%  undisturbed like you..." - Henri Cartan

