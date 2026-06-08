% Applies a propagator repeatedly by binary adaptive squaring. Syntax:
%
%              rho=multiprop(spin_system,P,rho,N)
%
% Parameters:
%
%     spin_system - Spinach spin system object
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
%       P^N explicitly. Propagator squares are cleaned up using
%       spin_system.tols.prop_zero.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=multiprop.m>

function rho=multiprop(spin_system,P,rho,N)

% Check consistency
grumble(spin_system,P,rho,N);

% Zero-step shortcut
if N==0, return; end

% Convert the step count into a bit-addressable integer
if ~isa(N,'uint64'), N=uint64(N); end

% Determine whether density-matrix action is needed
rho_is_matrix=~isvector(rho);

% Fast path for diagonal propagators
if isdiag(P)

    % Extract the diagonal of the propagator
    P_diag=diag(P);

    % Initialise the accumulated diagonal power
    P_power=ones(size(P_diag),'like',P_diag);

    % Raise the diagonal to the requested power
    while N>0

        % Accumulate the current binary power if present
        if bitand(N,uint64(1))>0
            P_power=P_diag.*P_power;
        end

        % Shift to the next binary digit
        N=bitshift(N,-1);

        % Square the diagonal only if higher powers remain
        if N>0, P_diag=clean_up(spin_system,P_diag.*P_diag,spin_system.tols.prop_zero); end

    end

    % Apply the diagonal power
    if rho_is_matrix
        rho=P_power.*rho.*P_power';
    else
        rho=P_power.*rho;
    end

    % Nothing else to do
    return;

end

% Run matrix-power accumulation for density matrices
if rho_is_matrix

    % Initialise the accumulated propagator
    P_total=[];

    % Process the binary expansion of the step count
    while N>0

        % Accumulate the current binary power if present
        if bitand(N,uint64(1))>0
            if isempty(P_total)
                P_total=P;
            else
                P_total=P*P_total;
            end
        end

        % Shift to the next binary digit
        N=bitshift(N,-1);

        % Square the propagator only if higher powers remain
        if N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_zero); end

    end

    % Apply the accumulated propagator
    rho=P_total*rho*P_total';

    % Nothing else to do
    return;

end

% Process the binary expansion of the step count
while N>0

    % Apply the current binary power if present
    if bitand(N,uint64(1))>0
        rho=P*rho;
    end

    % Shift to the next binary digit
    N=bitshift(N,-1);

    % Use two state-vector actions instead of the last matrix square
    if (N==1)&&(size(P,1)>1)
        rho=P*(P*rho);
        return;
    end

    % Square the propagator only if higher powers remain
    if N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_zero); end

end

end

% Consistency enforcement
function grumble(spin_system,P,rho,N)
if (~isstruct(spin_system))||(~isfield(spin_system,'tols'))||...
   (~isfield(spin_system.tols,'prop_zero'))
    error('spin_system.tols.prop_zero must exist.');
end
if (~isnumeric(spin_system.tols.prop_zero))||...
   (~isreal(spin_system.tols.prop_zero))||...
   (~isscalar(spin_system.tols.prop_zero))||...
   (spin_system.tols.prop_zero<0)
    error('spin_system.tols.prop_zero must be a non-negative real number.');
end
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

