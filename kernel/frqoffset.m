% Adds omega*Lz Zeeman frequency offsets to the Hamiltonian;
% this is useful in liquid state NMR experiments. Syntax:
%
%           H=frqoffset(spin_system,H,parameters)
%
% Parameters:
%
%                   H - Hamiltonian operator or commutati-
%                       on superoperator
%
%    parameters.spins - a cell array giving the spins that
%                       the offsetrs should be applied to,
%                       e.g. {'1H','13C'}
%
%   parameters.offset - a vector giving Zeeman offsets on
%                       each of the spins listed in pa-
%                       rameters.spins array.
%
% Outputs:
%
%                   H - Hamiltonian operator or commutati-
%                       on superoperator
%
% Note: offset transformation of this kind is an approximati-
%       on, use rotframe.m or intrep.m if a rigorous treat-
%       ment of second order effects is required.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=frqoffset.m>

function H=frqoffset(spin_system,H,parameters)

% Check consistency
grumble(spin_system,parameters);

% See if there are multiple channels on the same spin
[unique_spins,forward_index,backward_index]=unique(parameters.spins);

% Decide how to proceed
if numel(unique_spins)==numel(parameters.spins)
    
    % Simply apply the offsets
    for n=find(parameters.offset~=0)
        report(spin_system,['applying ' num2str(parameters.offset(n)) ...
                            ' Hz Zeeman frequency offset to ' parameters.spins{n} '...']);
        H=H+2*pi*parameters.offset(n)*operator(spin_system,'Lz',parameters.spins{n});
    end
    
else
    
    % Get unique offsets
    unique_offsets=parameters.offset(forward_index);
    
    % Check offsets on duplicate channels
    if ~all(unique_offsets(backward_index)==parameters.offset)
        error('offsets on channels referring to the same spin must be the same.');
    end
    
    % Apply the offsets
    for n=find(unique_offsets~=0)
        report(spin_system,['applying ' num2str(unique_offsets(n)) ...
                            ' Hz Zeeman frequency offset to ' unique_spins{n} '...']);
        H=H+2*pi*unique_offsets(n)*operator(spin_system,'Lz',unique_spins{n});
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,parameters)
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end
if ~isfield(parameters,'offset')
    error('parameters.offset variable must be present.');
end
if isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins must be specified alongside parameters.offset');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset must have the same number of elements as parameters.spins');
end
end

% [...] this flaw was identified by the brilliant German philosopher
% Friedrich Nietzsche who described it as "an inversion of morality"
% whereby the weak, the poor, the meek, the oppressed and the wretch-
% ed are virtuous and blessed by God whereas the strong, the wealthy,
% the noble and the powerful are the immoral and damned by the venge-
% ful almighty Yahweh for eternity. Nietzsche, with great insight and
% perception, stated that Christianity would be abandoned en masse in
% the twentieth century but that Westerners would still cling to this
% inversion of morality.
%
% Anders Behring Breivik

