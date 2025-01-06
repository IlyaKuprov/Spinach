% Scalar relaxation superoperator using Redfield theory. Syntax:
%
%            R=rlx_scalar(spin_system,H0,H1,tau_c_array)
%
% Parameters:
%
%     H0 - background Hamiltonian
%
%     H1 - the stochastically modulated interaction operator
%          multiplied by its root mean square modulation depth
%
%     tau_c_array - a cell array of the following format:
%
%                      {[weight_a,tau_a],[weight_b,tau_b],...}
%
%                   giving weights of the exponential components
%                   of the correlation function and the associa-
%                   ted correlation times, e.g. {[1.0,1e-12]}
%
% Outputs:
%
%     R  - relaxation superoperator, a negative definite matrix
%
% Note: if H1(t) has a non-zero time or ensemble average value,
%       that average must be subtracted out and placed into H0
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_scalar.m>

function R=rlx_scalar(spin_system,H0,H1,tau_c_array)

% Check consistency
grumble(H0,H1,tau_c_array);

% Get the ouput started
R=sparse(0);

% Loop over correlation function components
for n=1:numel(tau_c_array)
    
    % Extract component weight and correlation time
    weight=tau_c_array{n}(1); tau_c=tau_c_array{n}(2);
    
    % Ignore insignificant components
    if (weight~=0)&&(tau_c~=0)
    
        % Set the upper integration limit according to the accuracy goal
        upper_limit=2*tau_c*log(1/spin_system.tols.rlx_integration);

        % Remove inconsequential non-zeroes from H0
        H0=clean_up(spin_system,H0,1e-2/upper_limit);

    	% Take the integral using the auxiliary matrix exponential technique
        R=R-weight*H1*expmint(spin_system,H0,H1',H0+(1i/tau_c)*speye(size(H0)),upper_limit);
    
    end
    
end
    
end

% Consistency enforcement
function grumble(H0,H1,tau_c_array)
if (~isnumeric(H0))||(~isnumeric(H1))||(~ismatrix(H0))||...
   (~ismatrix(H1))||(~ishermitian(H0))||(~ishermitian(H1))
    error('H0 and H1 must be Hermitian square matrices.');
end
if ~all(size(H0)==size(H1))
    error('H0 and H1 must have the same dimension.');
end
if ~iscell(tau_c_array)
    error('tau_c_array must be a cell array of 2-element vectors.');
end
for n=1:numel(tau_c_array)
    if (~isnumeric(tau_c_array{n}))||(~isreal(tau_c_array{n}))||...
       (numel(tau_c_array{n})~=2)
        error('tau_c_array must be a cell array of 2-element vectors.');
    end
    if tau_c_array{n}(2)<0
        error('correlation times in tau_c_array must be non-negative.');
    end
end
end

% The dumbest rabbit is the one who thinks that, if he
% behaves himself, the wolves would not eat him.
%
% A Russian saying

