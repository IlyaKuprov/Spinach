% Computes integer propagator powers via an efficient powers-of-two
% based strategy. Syntax:
%
%              P=ppower(spin_system,P,N)
%
% Parameters:
%
%     spin_system - Spinach spin system object
%
%     P   - propagator matrix
%
%     N   - non-negative integer propagator power
%
% Outputs:
%
%     P   - propagator matrix raised to the power of N
%
% Note: the algorithm expands N into binary powers, squares P successively,
%       and multiplies only the active powers into the result. This avoids
%       explicit repeated multiplication. Propagator powers are cleaned up using
%       spin_system.tols.prop_chop.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ppower.m>

function P=ppower(spin_system,P,N)

% Check consistency
grumble(spin_system,P,N);

% Convert the step count into a bit-addressable integer
if ~isa(N,'uint64'), N=uint64(N); end

% Return the identity for the zero power
if N==0
    if issparse(P)
        P=speye(size(P,1));
    else
        P=eye(size(P,1),'like',P);
    end
    return
end

% First power shortcut
if N==1, return; end

% Initialise the propagator accumulator
if issparse(P)
    Q=speye(size(P,1));
else
    Q=eye(size(P,1),'like',P);
end

% Process the binary expansion of the step count
while N>0

    % Multiply the current binary power into the result if present
    if bitand(N,uint64(1))>0
        Q=clean_up(spin_system,Q*P,spin_system.tols.prop_chop);
    end

    % Shift to the next binary digit
    N=bitshift(N,-1);

    % Square the propagator only if higher powers remain
    if N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); end

end

% Return the accumulated propagator power
P=Q;

end

% Consistency enforcement
function grumble(spin_system,P,N)
if (~isstruct(spin_system))||(~isfield(spin_system,'tols'))||...
   (~isfield(spin_system.tols,'prop_chop'))
    error('spin_system.tols.prop_chop must exist.');
end
if (~isnumeric(spin_system.tols.prop_chop))||...
   (~isreal(spin_system.tols.prop_chop))||...
   (~isscalar(spin_system.tols.prop_chop))||...
   (spin_system.tols.prop_chop<0)
    error('spin_system.tols.prop_chop must be a non-negative real number.');
end
if (~isnumeric(P))||(~ismatrix(P))||(size(P,1)~=size(P,2))
    error('P must be a square matrix.');
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
if ~allfinite(P)
    error('P must be finite.');
end
end

% Your mind is software. Program it. Your body is a shell. Change it.
% Death is a disease. Cure it. Extinction is approaching. Fight it.
%
% Eclipse Phase


