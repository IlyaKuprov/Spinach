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
%     rho - state vector or state-vector stack in Liouville space or
%           wavefunction formalism, or a density matrix in Hilbert space
%           formalism
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
%       spin_system.tols.prop_chop.
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

% Determine propagation style from Spinach formalism
single_side=ismember(spin_system.bas.formalism,{'sphten-liouv',...
                                                'zeeman-liouv',...
                                                'zeeman-wavef'});

% Process the binary expansion of the step count
while N>0

    % Apply the current binary power if present
    if bitand(N,uint64(1))>0
        if single_side
            rho=P*rho;
        else
            rho=P*rho*P';
        end
    end

    % Shift to the next binary digit
    N=bitshift(N,-1);

    % Square the propagator only if higher powers remain
    if N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); end

end

end

% Consistency enforcement
function grumble(spin_system,P,rho,N)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
   (~isfield(spin_system.bas,'formalism'))
    error('spin_system.bas.formalism must exist.');
end
if (~ischar(spin_system.bas.formalism))||...
   (~ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv',...
                                         'sphten-liouv','zeeman-wavef'}))
    error('spin_system.bas.formalism is not recognised.');
end
if (~isfield(spin_system,'tols'))||...
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
if ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv',...
                                       'zeeman-wavef'})
    if size(P,1)~=size(rho,1)
        error('dimensions of P and rho must be consistent.');
    end
elseif size(rho,1)~=size(rho,2)
    error('rho must be a square matrix in Hilbert space formalism.');
elseif size(P,1)~=size(rho,1)
    error('dimensions of P and rho must be consistent.');
end
end

% Few forces can match the power of fanaticism. One
% that comes close is wounded pride.
%
% Frank Herbert, in
% the Dune series

