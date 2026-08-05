% Prints bosonic mode coupling summary for a Spinach system. This
% covers mode-mode exchange couplings, cross-Kerr couplings, spin-
% mode exchange couplings, longitudinal and radiation pressure
% couplings, and dispersive couplings. Syntax:
%
%               summary_mode_coup(spin_system,header)
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
% <https://spindynamics.org/wiki/index.php?title=summary_mode_coup.m>

function summary_mode_coup(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the summary table
report(spin_system,header);
report(spin_system,'====================================================');
report(spin_system,'Prt A   Prt B   Coupling type       Amplitude, Hz   ');
report(spin_system,'----------------------------------------------------');

% Print exchange couplings
[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.exchange));
for n=1:numel(rows)
    report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8)...
                        pad('exchange',20)...
                        num2str(spin_system.inter.modes.exchange{rows(n),cols(n)}/(2*pi),'%+10.4e')]);
end

% Print cross-Kerr couplings
[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.kerr));
for n=1:numel(rows)
    report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8)...
                        pad('cross-Kerr',20)...
                        num2str(spin_system.inter.modes.kerr{rows(n),cols(n)}/(2*pi),'%+10.4e')]);
end

% Print longitudinal and radiation pressure couplings
[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.longitudinal));
for n=1:numel(rows)
    if all(ismember(spin_system.comp.types([rows(n) cols(n)]),{'C','V','T'}))
        coup_name='radiation pressure';
    else
        coup_name='longitudinal';
    end
    report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8)...
                        pad(coup_name,20)...
                        num2str(spin_system.inter.modes.longitudinal{rows(n),cols(n)}/(2*pi),'%+10.4e')]);
end

% Print dispersive couplings
[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.dispersive));
for n=1:numel(rows)
    report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8)...
                        pad('dispersive',20)...
                        num2str(spin_system.inter.modes.dispersive{rows(n),cols(n)}/(2*pi),'%+10.4e')]);
end
report(spin_system,'====================================================');

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


