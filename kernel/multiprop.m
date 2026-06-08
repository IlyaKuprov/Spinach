% Applies a propagator repeatedly by binary adaptive squaring. Syntax:
%
%                     rho=multiprop(P,rho,N)
%
% Parameters:
%
%     P   - propagator matrix
%
%     rho - state vector in Liouville space or wavefunction formalism,
%           or a density matrix in Hilbert space formalism
%
%     N   - number of times to apply the propagator
%
% Outputs:
%
%     rho - state vector or density matrix after N applications of P
%
% Note: the algorithm expands N into binary powers, squares P successively,
%       and applies only the active powers to rho. This avoids constructing
%       P^N explicitly.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=multiprop.m>

function rho=multiprop(P,rho,N)

% Check consistency
grumble(P,rho,N);

% Zero-step shortcut
if N==0, return; end

% Convert the step count into a bit-addressable integer
if ~isa(N,'uint64'), N=uint64(N); end

% Determine whether density-matrix action is needed
rho_is_matrix=~isvector(rho);

% Process the binary expansion of the step count
while N>0

    % Apply the current binary power if present
    if bitand(N,uint64(1))>0
        if rho_is_matrix
            rho=P*rho*P';
        else
            rho=P*rho;
        end
    end

    % Shift to the next binary digit
    N=bitshift(N,-1);

    % Square the propagator only if higher powers remain
    if N>0, P=P*P; end

end

end

% Consistency enforcement
function grumble(P,rho,N)
if (~isnumeric(P))||(~ismatrix(P))||(size(P,1)~=size(P,2))
    error('P must be a square matrix.');
end
if (~isnumeric(rho))||(~ismatrix(rho))
    error('rho must be a numeric matrix.');
end
if (~isnumeric(N))||(~isscalar(N))||(~isreal(N))
    error('N must be a non-negative real integer.');
end
if isinteger(N)
    if N<0
        error('N must be a non-negative real integer.');
    end
else
    if (~isfinite(N))||(N<0)||(mod(N,1)~=0)||(N>flintmax)
        error('N must be a non-negative real integer not exceeding flintmax.');
    end
end
if (~allfinite(P))||(~allfinite(rho))
    error('P and rho must be finite.');
end
if isvector(rho)
    if ~iscolumn(rho)
        error('rho must be a column vector or a square matrix.');
    end
    if size(P,1)~=size(rho,1)
        error('dimensions of P and rho must be consistent.');
    end
elseif size(rho,1)~=size(rho,2)
    error('rho must be a column vector or a square matrix.');
elseif size(P,1)~=size(rho,1)
    error('dimensions of P and rho must be consistent.');
end
end

% Your mind is software. Program it. Your body is a shell. Change it.
% Death is a disease. Cure it. Extinction is approaching. Fight it.
%
% Eclipse Phase

