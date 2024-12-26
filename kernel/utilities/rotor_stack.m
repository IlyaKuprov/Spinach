% Returns a rotor stack of Liouvillians or Hamiltonians. The stack is 
% needed for the traditional style calculation of MAS dynamics. Syntax:
%
%            L=rotor_stack(spin_system,parameters,assumptions)
%
% Parameters:
%
%   parameters.axis     - spinning axis, given as a normalized
%                         3-element vector
%
%   parameters.offset   - a cell array giving transmitter off-
%                         sets in Hz on each of the spins listed
%                         in parameters.spins array
%
%   parameters.spins    - a cell array giving the spins that 
%                         the offsets refer to, e.g. {'1H','13C'}
%
%   parameters.max_rank - maximum harmonic rank to retain in
%                         the solution (increase till conver-
%                         gence is achieved, approximately
%                         equal to the number of spinning si-
%                         debands in the spectrum)
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
%   parameters.orientation - the orientation of the spin system
%                            at rotor phase zero, a vector of 
%                            three Euler angles in radians.
%
%   parameters.masframe    - the frame in which the rotations 
%                            are applied. The possibilities are:
%
%            'magnet' - the initial orientation in the lab frame
%                       (three-angle powder grids will be required)
%
%            'rotor'  - the initial orientation in the rotor frame
%                       (two-angle powder grids will be required)
%
%   assumptions - assumption set to be used in generating the
%                 Hamiltonian, see assume.m
%
% Outputs:
%
%   L   - a cell array of Hamiltonian or Liouvillian matrices,
%         one for each tick of the rotor.
%
%   rotor_phases - rotor phases at each tick, radians
%
% Note: relaxation and chemical kinetics are not included.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rotor_stack.m>

function [L,rotor_phases]=rotor_stack(spin_system,parameters,assumptions)

% Check consistency
grumble(spin_system,parameters,assumptions);

% Get the Hamiltonian
[H,Q]=hamiltonian(assume(spin_system,assumptions));

% Apply offsets
H=frqoffset(spin_system,H,parameters);

% Get rotor axis orientation
[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1),...
                                   parameters.axis(2),...
                                   parameters.axis(3)); 
rotor_theta=pi/2-rotor_theta;

% Compute rotor angles
rotor_phases=fourdif(2*parameters.max_rank+1,1);

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Preallocate Liouvillian blocks
L=cell(2*parameters.max_rank+1,1);

% Build Liouvillian blocks
parfor n=1:(2*parameters.max_rank+1) %#ok<*PFBNS>
    
    % Get the block started
    L{n}=H;
    
    % Loop over spherical ranks
    for r=1:numel(Q)
        
        % Compute rotor axis tilt
        D_rot2lab=wigner(r,+rotor_phi,+rotor_theta,0);
        D_lab2rot=wigner(r,0,-rotor_theta,-rotor_phi);
        
        % Compute the initial crystallite orientation
        D_initial=wigner(r,parameters.orientation(1),...
                           parameters.orientation(2),...
                           parameters.orientation(3));
    
         % Compute rotor rotation
        D_rotor=wigner(r,0,0,rotor_phases(n));
        
        % Decide the reference frame
        if strcmp(parameters.masframe,'magnet')
            
            % Compose rotations
            D=D_rot2lab*D_rotor*D_lab2rot*D_initial;
            
        elseif strcmp(parameters.masframe,'rotor')
            
            % Compose rotations
            D=D_rot2lab*D_rotor*D_initial;
            
        else
            
            % Complain and bomb out
            D=0; error('unknown MAS frame.'); %#ok<NASGU>
            
        end
        
        % Build the block
        for k=1:(2*r+1)
            for m=1:(2*r+1)
                L{n}=L{n}+D(k,m)*Q{r}{k,m};
            end
        end
        
    end
    
    % Apply interaction representations
    for k=1:numel(parameters.rframes)
        L{n}=rotframe(spin_system,C{k},(L{n}+L{n}')/2,parameters.rframes{k}{1},...
                                                      parameters.rframes{k}{2});
    end
    
    % Clean up the result
    L{n}=clean_up(spin_system,L{n},spin_system.tols.liouv_zero);
    
end

end

% Consistency enforcement
function grumble(spin_system,parameters,assumptions)

% Spinning axis
if ~isfield(parameters,'axis')
    error('parameters.axis subfield must be present.');
elseif (~isnumeric(parameters.axis))||(~isreal(parameters.axis))||...
       (~isrow(parameters.axis))||(numel(parameters.axis)~=3)
    error('parameters.axis must be a row vector of three real numbers.');
end

% Offsets
if ~isfield(parameters,'offset')
    error('offsets must be specified in parameters.offset field.');
elseif isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins variable must be specified alongside parameters.offset variable.');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset variable must have the same number of elements as parameters.spins.');
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

% Rotor rank
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank subfield must be present.');
elseif (~isnumeric(parameters.max_rank))||(~isreal(parameters.max_rank))||...
       (mod(parameters.max_rank,1)~=0)||(parameters.max_rank<0)
    error('parameters.max_rank must be a positive real integer.');
end

% Initial orientation
if ~isfield(parameters,'orientation')
    error('initial orientation must be specified in parameters.orientation variable.');
elseif ~all(size(parameters.orientation)==[1 3])
    error('parameters.orientation variable must be a row vector with three numbers.');
elseif ~isnumeric(parameters.orientation)
    error('elements of parameters.orientation vector must be real numbers.');
end
if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
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

% MAS frame
if ~isfield(parameters,'masframe')
    error('powder averaging frame must be specified in parameters.masframe field.');
elseif ~ischar(parameters.masframe)
    error('parameters.masframe must be a character string.');
elseif ~ismember(parameters.masframe,{'rotor','magnet'})
    error('the value of parameters.masframe is not recognised.');
end
    
end

% Never complain and never explain. 
%
% Benjamin Disraeli

