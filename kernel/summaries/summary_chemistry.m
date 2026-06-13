% Prints chemical subsystem and exchange summary for a Spinach system. Syntax:
%
%                 summary_chemistry(spin_system)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object
%
% Outputs:
%
%    this function prints to the console or to the user-specified
%    output via report.m function
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=summary_chemistry.m>

function summary_chemistry(spin_system)

% Check consistency
grumble(spin_system);

% Report multiple chemical subsystems
if numel(spin_system.chem.parts)>1

    % Report spin system partitioning
    for n=1:numel(spin_system.chem.parts)
        report(spin_system,['chemical subsystem ' num2str(n) ' contains spins: ' num2str(spin_system.chem.parts{n})]);
    end

    % Report first-order reaction rates
    if isfield(spin_system.chem,'rates')
        report(spin_system,'inter-subsystem reaction rates:');
        report(spin_system,'===============================');
        report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
        report(spin_system,'-------------------------------');
        [rows,cols,vals]=find(spin_system.chem.rates);
        for n=1:length(vals)
            report(spin_system,[' ' strjust([num2str(rows(n)) blanks(3-length(num2str(rows(n))))],'left') '       '...
                                    strjust([num2str(cols(n)) blanks(3-length(num2str(cols(n))))],'left') '      '...
                                             num2str(vals(n),'%+0.3e')]);
        end
        report(spin_system,'===============================');
    end

end

% Report flux rates if specified
if isfield(spin_system.chem,'flux_rate')
    [rows,cols,vals]=find(spin_system.chem.flux_rate);
    if numel(vals)>0
        report(spin_system,'point-to-point flux rates:');
        report(spin_system,'===============================');
        report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
        report(spin_system,'-------------------------------');
        for n=1:length(vals)
            report(spin_system,[' ' strjust([num2str(rows(n)) blanks(3-length(num2str(rows(n))))],'left') '       '...
                                    strjust([num2str(cols(n)) blanks(3-length(num2str(cols(n))))],'left') '      '...
                                             num2str(vals(n),'%+0.3e')]);
        end
    end
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

% Always code as if the guy who ends up maintaining 
% your code will be a violent psychopath who knows 
% where you live.
%
% Martin Golding

