% Prints bosonic mode parameter summary for a Spinach system. Syntax:
%
%                 summary_modes(spin_system,header)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object
%
%    header       - a string of text to precede the summary
%
% Outputs:
%
%    this function prints to the console or to the user-specified
%    output via report.m function
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=summary_modes.m>

function summary_modes(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the summary table
report(spin_system,header);
report(spin_system,'==========================================================================================');
report(spin_system,'#    Mode  Type      Levels  Freq, Hz      Anharm, Hz    Damping, Hz   Dephasing, Hz     ');
report(spin_system,'------------------------------------------------------------------------------------------');
for n=1:spin_system.comp.nspins
    if ismember(spin_system.comp.types{n},{'C','V','T'})

        % Translate the type letter into a word
        switch spin_system.comp.types{n}
            case 'C', type_word='cavity';
            case 'V', type_word='phonon';
            case 'T', type_word='transmon';
        end

        % Do the printing in Hz
        report(spin_system,[pad(num2str(n),5)...
                            pad(spin_system.comp.isotopes{n},6)...
                            pad(type_word,10)...
                            pad(num2str(spin_system.comp.mults(n)),8)...
                            pad(num2str(spin_system.inter.modes.frqs(n)/(2*pi),'%+10.4e'),14)...
                            pad(num2str(spin_system.inter.modes.anharms(n)/(2*pi),'%+10.4e'),14)...
                            pad(num2str(spin_system.inter.modes.damp(n)/(2*pi),'%+10.4e'),14)...
                            pad(num2str(spin_system.inter.modes.dephase(n)/(2*pi),'%+10.4e'),14)]);

    end
end
report(spin_system,'==========================================================================================');

end

% Consistency enforcement
function grumble(spin_system,header)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if ~ischar(header)
    error('header must be a character string.');
end
end

