% Converts directional cosine matrix into Euler angles, ZYZ active
% convention (rotating the object rather than the axes). Syntax:
%
%                [alpha,beta,gamma]=dcm2euler(dcm)
%
%                               OR
%
%                      angles=dcm2euler(dcm)
%
% Parameters:
%
%     dcm                - directional cosine matrix
%
% Outputs:
%
%     alpha, beta, gamma - Euler angles in ZYZ active con-
%                          vention, radians
%
%     angles             - a row vector of Euler angles in
%                          ZYZ active convention, ordered 
%                          as alpha, beta, gamma, in radians
%
% Note: the problem of recovering Euler angles from a DCM is, in 
%       general, ill-posed. This function is a product of consi-
%       derable work, it has passed rigorous testing: it either
%       returns a correct answer or gives an informative error.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dcm2euler.m>

function [arg1,arg2,arg3]=dcm2euler(dcm)

% Check consistency
grumble(dcm);

% Get the beta angle out and wrap it into [0,pi]
beta=mod(acos(dcm(3,3)),pi);

% Do a brute force surface scan with respect to alpha and gamma
alphas=pi*linspace(0.05,1.95,20);
gammas=pi*linspace(0.05,1.95,20);
n_min=1; k_min=1; err_min=1;
for n=1:20
    for k=1:20
        err_current=norm(euler2dcm(alphas(n),beta,gammas(k))-dcm,1);
        if err_current<err_min
            n_min=n; k_min=k; err_min=err_current;
        end
    end
end
alpha=alphas(n_min); gamma=gammas(k_min);

% Run the optimisation on alpha and gamma
options=optimset('Display','off','LargeScale','off','TolX',1e-12,'TolFun',1e-12);
answer=fminunc(@(angles)norm(euler2dcm(angles(1),beta,angles(2))-dcm,'fro')^2,[alpha gamma],options);
alpha=answer(1); gamma=answer(2);

% Wrap both angles into [0,2*pi]
alpha=mod(alpha,2*pi); gamma=mod(gamma,2*pi);

% Make sure the result is good enough and bomb out if not
if norm(dcm-euler2dcm(alpha,beta,gamma),1)>1e-3
    disp(dcm); disp(euler2dcm(alpha,beta,gamma));
    error('DCM to Euler conversion failed.');
end

% Adapt to the output style
if nargout==1||nargout==0
    arg1=[alpha beta gamma];
elseif nargout==3
    arg1=alpha; arg2=beta; arg3=gamma;
else
    error('incorrect number of output arguments.');
end

end

% Consistency enforcement
function grumble(dcm)
if (~isnumeric(dcm))||(~isreal(dcm))||(~all(size(dcm)==[3 3]))
    error('DCM must be a real 3x3 matrix.');
end
if norm(dcm'*dcm-eye(3),1)>1e-6
    warning('DCM is not orthogonal to 1e-6 tolerance - conversion accuracy not guaranteed.');
end
if norm(dcm'*dcm-eye(3),1)>1e-2
    error('DCM is not orthogonal to 1e-2 tolerance, cannot proceed with conversion.');
end
if abs(det(dcm)-1)>1e-6
    warning('DCM determinant is not unit to 1e-6 tolerance - conversion accuracy not guaranteed.');
end
if abs(det(dcm)-1)>1e-2
    error('DCM determinant is not unit to 1e-2 tolerance, cannot proceed with conversion.');
end
end

% I would give the greatest sunset in the world for one sight of New York's
% skyline. Particularly when one can't see the details. Just the shapes. The
% shapes and the thought that made them. The sky over New York and the will
% of man made visible. What other religion do we need? And then people tell
% me about pilgrimages to some dank pesthole in a jungle where they go to do
% homage to a crumbling temple, to a leering stone monster with a pot belly,
% created by some leprous savage. Is it beauty and genius they want to see?
% Do they seek a sense of the sublime? Let them come to New York, stand on
% the shore of the Hudson, look and kneel. When I see the city from my win-
% dow - no, I don't feel how small I am - but I feel that if a war came to
% threaten this, I would like to throw myself into space, over the city, and
% protect these buildings with my body.
%
% Ayn Rand, "The Fountainhead"

