% Single-crystal interface to pulse sequences. Generates a Liouvillian
% superoperator and passes it on to the pulse sequence function, which
% should be supplied as a handle. Syntax:
%
%   answer=crystal(spin_system,pulse_sequence,parameters,assumptions)
%
% where pulse sequence is a function handle to one of the pulse sequences
% located in the experiments directory, assumptions is a string that would
% be passed to assume.m when the Hamiltonian is built and parameters is a
% structure with the following subfields:
%
%   parameters.spins    - a cell array giving the spins that 
%                         the pulse sequence involves, e.g. 
%                         {'1H','13C'}
%
%   parameters.offset   - a cell array giving transmitter off-
%                         sets in Hz on each of the spins listed
%                         in parameters.spins array
%
%   parameters.orientation - a row vector of the three Euler angles
%                           (in radians) giving the orientation of
%                            the system relative to the input orien-
%                            tation.
%
%   parameters.rframes  - rotating frame specification, e.g.
%                         {{'13C',2},{'14N,3}} requests second
%                         order rotating frame transformation
%                         with respect to carbon-13 and third
%                         order rotating frame transformation
%                         with respect to nitrogen-14. When
%                         this option is used, the assumptions
%                         on the respective spins should be
%                         laboratory frame.
%
%  parameters.needs   - a cell array of strings specifying additional
%                       information required by the sequence:
%
%                       'zeeman_op' - Zeeman part of the Hamiltonian
%                       in the laboratory frame, to be placed into
%                       parameters.hzeeman and sent to pulse sequence
%
%                       'aniso_eq' - thermal equilibrium is recomputed 
%                       using the full anisotropic Hamiltonian at the
%                       current orientation, and sent to the pulse 
%                       sequence in parameters.rho0 subfield
%
%  parameters.*       - additional subfields may be required by your
%                       pulse sequence - check its documentation page 
%
% The parameters structure is passed to the pulse sequence with the follo-
% wing additional parameters set:
%
%   parameters.spc_dim  - matrix dimension for the spatial
%                         dynamics subspace (1 in this case)
%
%   parameters.spn_dim  - matrix dimension for the spin 
%                         dynamics subspace
%
% This function returns whatever it is that the pulse sequence returns.
%
% Note: arbitrary order rotating frame transformation is supported, inc-
%       luding infinite order. See the header of rotframe.m for further
%       information.
%
% i.kuprov@soton.ac.nn.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=crystal.m>

function answer=crystal(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Set the assumptions
spin_system=assume(spin_system,assumptions);

% Get the rotational expansion
[I,Q]=hamiltonian(spin_system);

% Get the lab frame Zeeman operator if needed
if ismember('zeeman_op',parameters.needs)
    report(spin_system,'building the lab frame Zeeman operator...');
    [ZI,ZQ]=hamiltonian(assume(spin_system,'labframe','zeeman'));
end

% Compute the thermal equilibrium at the current orientation
if ismember('aniso_eq',parameters.needs)
    [HL,QL]=hamiltonian(assume(spin_system,'labframe'),'left');
    parameters.rho0=equilibrium(spin_system,HL,QL,parameters.orientation);
end

% Apply the offsets
I=frqoffset(spin_system,I,parameters);

% Build the Hamiltonian and tidy up
H=I+orientation(Q,parameters.orientation);
H=(H+H')/2; clear('I','Q');

% Get the lab frame Zeeman operator at the current orientation
if ismember('zeeman_op',parameters.needs)
    Z=ZI+orientation(ZQ,parameters.orientation);
    parameters.hzeeman=(Z+Z')/2; clear('ZI','ZQ');
end

% Apply rotating frames
for k=1:numel(parameters.rframes)
    
    % Get the carrier operator
    C=carrier(spin_system,parameters.rframes{k}{1});
    
    % Compute the rotating frame transformation
    H=rotframe(spin_system,C,H,parameters.rframes{k}{1},...
                               parameters.rframes{k}{2});
                           
end

% Build relaxation at the orientation specified
R=relaxation(spin_system,parameters.orientation);

% Build kinetics
K=kinetics(spin_system);

% Get problem dimensions
parameters.spc_dim=1; parameters.spn_dim=size(H,1);

% Report to the user
report(spin_system,'running the pulse sequence...');

% Call the pulse sequence
answer=pulse_sequence(spin_system,parameters,H,R,K);

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'decouple')
    report(spin_system,'parameters.decouple field not set, assuming no decoupling.');
    parameters.decouple={};
end
if ~isfield(parameters,'rframes')
    report(spin_system,'parameters.rframes field not set, assuming no additional rotating frames.');
    parameters.rframes={};
end
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'verbose')
    report(spin_system,'parameters.verbose field not set, silencing array operations.');
    parameters.verbose=0;
end
if ~isfield(parameters,'needs')
    report(spin_system,'parameters.needs field not set, assumpting empty.');
    parameters.needs={};
end
end

% Consistency checking
function grumble(spin_system,pulse_sequence,parameters,assumptions)

% Orientation
if ~isfield(parameters,'orientation')
    error('system orientation must be specified in parameters.orientation variable.');
elseif ~all(size(parameters.orientation)==[1 3])
    error('parameters.orientation variable must be a row vector with three numbers.');
elseif ~isnumeric(parameters.orientation)
    error('elements of parameters.orientation vector must be real numbers.');
end

% Pulse sequence
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end

% Assumptions
if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
end

% Active spins
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end

% Offsets
if isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins variable must be specified alongside parameters.offset variable.');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset variable must have the same number of elements as parameters.spins.');
end

% Rotating frame transformations
if ~isfield(parameters,'rframes')
    error('parameters.rframes variable must be specified.');
elseif ~iscell(parameters.rframes)
    error('parameters.rframes must be a cell array.');
end
for n=1:numel(parameters.rframes)
    if ~iscell(parameters.rframes{n})
        error('elements of parameters.rframes must be cell arrays.');
    end
    if numel(parameters.rframes{n})~=2
        error('elements of parameters.rframes must have exactly two sub-elements each.');
    end
    if ~ischar(parameters.rframes{n}{1})
        error('the first part of each element of parameters.rframes must be a character string.');
    end
    if ~ismember(parameters.rframes{n}{1},spin_system.comp.isotopes)
        error('parameters.rframes refers to a spin that is not present in the system.');
    end
    if ~isnumeric(parameters.rframes{n}{2})
        error('the second part of each element of parameters.rframes must be a number.');
    end
end

end

% Why do they always teach us that it's easy and evil to do what we 
% want and that we need discipline to restrain ourselves? It's the
% hardest thing in the world - to do what we want. And it takes the
% greatest kind of courage. I mean, what we really want.
%
% Ayn Rand, "The Fountainhead"

