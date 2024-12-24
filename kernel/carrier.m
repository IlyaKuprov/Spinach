% Returns the "carrier" Hamiltonian - the part of the Zeeman interaction
% Hamiltonian that corresponds to all particles having the Zeeman frequ-
% ency prescribed by their isotropic free-particle magnetogyric ratio and
% Z axis magnet field specified by the user. This Hamiltonian is used in
% rotating frame transforms and average Hamiltonian theories. Syntax:
%
%              H=carrier(spin_system,spins,operator_type)
%
% Parameters:
%
%     spins - a string specifying the isotope, e.g. '1H';
%             to select all spins, use 'all'.
%
%     in Liouville space, operator_type can be set to
%                      
%           'left' - produces left side product superoperator
%
%          'right' - produces right side product superoperator
%
%           'comm' - produces commutation superoperator (default)
%
%          'acomm' - produces anticommutation superoperator
%
%     in Hilbert space this parameter is ignored.
%
% Outputs:
%
%     H - a Hamiltonian (Hilbert space) or its superoperator 
%         of the specified type (Liouville space).
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=carrier.m>

function H=carrier(spin_system,spins,operator_type)

% Set the default for the type
if ~exist('operator_type','var'), operator_type='comm'; end

% Check consistency
grumble(spin_system,spins);

% Preallocate the answer
H=mprealloc(spin_system,1);

% Find the spins
if strcmp(spins,'all')
    spin_index=1:spin_system.comp.nspins;
else
    spin_index=find(strcmp(spins,spin_system.comp.isotopes));
end

% Compute the answer
for n=1:numel(spin_index)
    if abs(spin_system.inter.basefrqs(spin_index(n)))>0
        H=H+spin_system.inter.basefrqs(spin_index(n))*...
            operator(spin_system,{'Lz'},{spin_index(n)},operator_type);
    end
end

% Clean up the answer
H=clean_up(spin_system,(H+H')/2,spin_system.tols.liouv_zero);

end

% Consistency enforcement
function grumble(spin_system,spins)
if ~ischar(spins)
    error('spins argument must be a character array');
end
if (~strcmp(spins,'all'))&&(~ismember(spins,spin_system.comp.isotopes))
    error('no such spins in the system.');
end
end

% Atelophobia - (n.) fear of imperfection

