% Generates a Lindblad superoperator from user-specified left-side and
% right-side product superoperators and calibrates it using the experi-
% mental relaxation rate of a user-specified state. Syntax:
%
%              R=lindbladian(A_left,A_right,rho,rlx_rate)
%
% Parameters:
%
%     A_left   - left side product superoperator of the 
%                interaction that is causing relaxation
%                (see operator.m and hamiltonian.m)
%
%     A_right  - right side product superoperator of the
%                same interaction
%
%     rho      - the state vector whose relaxation rate
%                is known from the experiment 
%
%     rlx_rate - experimental relaxation rate of rho
%
% Outputs:
%
%     R        - Lindblad relaxation superoperator indu-
%                ced by the interaction A, such that
%                <rho|R|rho>/norm(rho,2)^2 = -rlx_rate
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lindbladian.m>

function R=lindbladian(A_left,A_right,rho,rlx_rate)

% Check consistency
grumble(A_left,A_right,rho,rlx_rate);

% Generate a Lindbladian
R=A_left*A_right'-(A_left'*A_left+A_right*A_right')/2;

% Check for silly inputs
if abs(rho'*R*rho)<1e-10
    error('the operator supplied does not appear to relax the given state.');
end

% Calibrate the Lindbladian
rho=rho/norm(rho,2); R=-rlx_rate*R/(rho'*R*rho);

end

% Consistency enforcement
function grumble(A_left,A_right,rho,rlx_rate)
if (~isnumeric(A_left))||(~isnumeric(A_right))||...
   (~isnumeric(rho))||(~isnumeric(rlx_rate))
    error('all inputs must be numeric.');
end
if (~isnumeric(rlx_rate))||(~isscalar(rlx_rate))||...
   (~isreal(rlx_rate))
    error('rlx_rate must be a non-negative real number.');
end
end

% Morality, it could be argued, represents the way that people
% would like the world to work, wheareas economics represents
% how it actually does work.
%
% Steven D. Levitt, "Freakonomics"

