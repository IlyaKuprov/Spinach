% Spin-boson device interface to pulse sequences. Generates the evolution
% generators for a device containing spins and bosonic modes at a fixed
% orientation of the spin subsystem, and passes them to the pulse sequence
% function, which should be supplied as a handle. Syntax:
%
%   answer=device(spin_system,pulse_sequence,parameters,assumptions)
%
% Parameters:
%
%  pulse_sequence         - pulse sequence function handle. See the
%                           experiments directory for the list of
%                           pulse sequences that ship with Spinach.
%
%  parameters.spins       - a cell array giving the spin species that
%                           the pulse sequence works on, in the order
%                           of channels, e.g. {'E'}; may be omitted
%                           when no spin channels are needed
%
%  parameters.offset      - transmitter offsets in Hz on each of the
%                           spin species listed in parameters.spins
%
%  parameters.mode_offset - detuning offsets in Hz, one for each
%                           bosonic mode in the order of declaration;
%                           the transmitter sign convention of para-
%                           meters.offset applies: each offset enters
%                           the Hamiltonian as minus the offset times
%                           the number operator of its mode, so that a
%                           positive offset lowers the mode frequency
%
%  parameters.decouple    - a cell array of spin species to be wiped
%                           from the evolution generators and from the
%                           initial state, e.g. {'1H'}; the default is
%                           an empty cell array, meaning no decoupling
%
%  parameters.orientation - Euler angles (ZYZ active convention, ra-
%                           dians) giving the orientation of the spin
%                           subsystem; bosonic terms are not affected;
%                           the default is [0 0 0]
%
%  parameters.rframes     - numerical rotating frame specification for
%                           spin species, e.g. {{'E',2}}, as described
%                           in the header of rotframe.m
%
%  parameters.needs       - a cell array of strings specifying additional
%                           information required by the sequence:
%
%                           'rho_eq' - thermal equilibrium state at the
%                           system temperature, with the Bose-Einstein
%                           populations of the bosonic modes included;
%                           this is placed into parameters.rho0
%
%  parameters.*           - additional subfields may be required by
%                           the pulse sequence - check its docs
%
%  assumptions            - 'labframe', 'cavity', or 'spin-phonon';
%                           see the header of assume.m for details
%
% Outputs:
%
%  answer - whatever it is that the pulse sequence returns.
%
% Note: the system must contain at least one bosonic mode; pure spin
%       systems belong with liquid.m, crystal.m, and powder.m contexts.
%
% Note: dissipative bosonic modes require a Liouville space formalism;
%       coherent simulations may also use zeeman-hilb.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=device.m>

function answer=device(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner');

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Set assumptions
spin_system=assume(spin_system,assumptions);

% Get the Hamiltonian at the spin subsystem orientation
[I,Q]=hamiltonian(spin_system);
H=I+orientation(Q,parameters.orientation);
H=(H+H')/2;

% Get relaxation and kinetics superoperators
R=relaxation(spin_system,parameters.orientation);
K=kinetics(spin_system);

% Get the thermal equilibrium if needed
if ismember('rho_eq',parameters.needs)
    report(spin_system,'building the thermal equilibrium state...');
    [HL,QL]=hamiltonian(assume(spin_system,'labframe'),'left');
    parameters.rho0=equilibrium(spin_system,HL,QL,parameters.orientation);
end

% Process spin channel offsets
if ~isempty(parameters.spins)
    H=frqoffset(spin_system,H,parameters);
end

% Process mode detuning offsets
mode_list=find(ismember(spin_system.comp.types,{'C','V','T'}));
for n=1:numel(mode_list)
    if abs(parameters.mode_offset(n))>0
        report(spin_system,['detuning offset for bosonic mode ' num2str(mode_list(n)) ...
                            ': ' num2str(parameters.mode_offset(n)) ' Hz']);
        H=H-2*pi*parameters.mode_offset(n)*operator(spin_system,{'N'},{mode_list(n)});
    end
end

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Apply rotating frames
for k=1:numel(parameters.rframes)
    H=rotframe(spin_system,C{k},H,parameters.rframes{k}{1},parameters.rframes{k}{2});
end

% Wipe the states of the decoupled spins
H=decouple(spin_system,H,[],parameters.decouple);
R=decouple(spin_system,R,[],parameters.decouple);
K=decouple(spin_system,K,[],parameters.decouple);
if ismember('rho_eq',parameters.needs)
    [~,parameters.rho0]=decouple(spin_system,[],parameters.rho0,parameters.decouple);
end

% Get problem dimensions
parameters.spc_dim=1; parameters.spn_dim=size(H,1);

% Report to the user
report(spin_system,'running the pulse sequence...');

% Call the pulse sequence
answer=pulse_sequence(spin_system,parameters,H,R,K);

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'spins')
    report(spin_system,'parameters.spins field not set, assuming no spin channels.');
    parameters.spins={};
end
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'mode_offset')
    report(spin_system,'parameters.mode_offset field not set, assuming zero mode offsets.');
    parameters.mode_offset=zeros(1,nnz(ismember(spin_system.comp.types,{'C','V','T'})));
end
if ~isfield(parameters,'orientation')
    report(spin_system,'parameters.orientation field not set, using the input orientation.');
    parameters.orientation=[0 0 0];
end
if ~isfield(parameters,'rframes')
    report(spin_system,'parameters.rframes field not set, assuming no additional rotating frame transformations.');
    parameters.rframes={};
end
if ~isfield(parameters,'needs')
    report(spin_system,'parameters.needs field not set, assuming empty.');
    parameters.needs={};
end
if ~isfield(parameters,'decouple')
    report(spin_system,'parameters.decouple field not set, assuming no decoupling.');
    parameters.decouple={};
end
end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters,assumptions)

if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end

if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
end
if ~ismember(assumptions,{'labframe','cavity','spin-phonon'})
    error('assumptions must be labframe, cavity, or spin-phonon in this context.');
end

if ~any(ismember({'C','V','T'},spin_system.comp.types),'all')
    error('no bosonic modes in the system, use liquid, crystal, or powder contexts.');
end

if isfield(spin_system.inter,'modes')&&...
   (~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'}))&&...
   any([spin_system.inter.modes.damp spin_system.inter.modes.dephase]>0)
    error('dissipative bosonic modes require a Liouville space formalism.');
end

if ~iscell(parameters.spins)
    error('parameters.spins must be a cell array of character strings.');
end
for n=1:numel(parameters.spins)
    if ~ischar(parameters.spins{n})
        error('all elements of parameters.spins must be character strings.');
    end
    if ~ismember(parameters.spins{n},spin_system.comp.isotopes)
        error('parameters.spins refers to a particle that is not present in the system.');
    end
    particles=strcmp(spin_system.comp.isotopes,parameters.spins{n});
    if ~all(strcmp(spin_system.comp.types(particles),'S'))
        error('parameters.spins may only contain spin species.');
    end
end

if ~isnumeric(parameters.offset)
    error('parameters.offset must be an array of real numbers.');
end
if numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset must have the same number of elements as parameters.spins.');
end

n_modes=nnz(ismember(spin_system.comp.types,{'C','V','T'}));
if (~isnumeric(parameters.mode_offset))||(~isreal(parameters.mode_offset))||...
   any(~isfinite(parameters.mode_offset),'all')
    error('parameters.mode_offset must be an array of finite real numbers.');
end
if numel(parameters.mode_offset)~=n_modes
    error('parameters.mode_offset must have one element per bosonic mode.');
end

if (~isnumeric(parameters.orientation))||(~isreal(parameters.orientation))||...
   (numel(parameters.orientation)~=3)
    error('parameters.orientation must be a real three-element vector of Euler angles.');
end

if ~iscell(parameters.rframes)
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
        error('parameters.rframes refers to a particle that is not present in the system.');
    end
    particles=strcmp(spin_system.comp.isotopes,parameters.rframes{n}{1});
    if ~all(strcmp(spin_system.comp.types(particles),'S'))
        error('parameters.rframes may only refer to spin species.');
    end
    if ~isnumeric(parameters.rframes{n}{2})
        error('the second part of each element of parameters.rframes must be a number.');
    end
end

if (~iscell(parameters.needs))||any(~cellfun(@ischar,parameters.needs(:)))
    error('parameters.needs must be a cell array of character strings.');
end
if any(~ismember(parameters.needs(:),{'rho_eq'}))
    error('parameters.needs may only contain ''rho_eq'' in this context.');
end

end

% There is nothing so practical as a good theory.
% Kurt Lewin


