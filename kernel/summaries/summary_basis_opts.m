% Prints basis-set option summary for a Spinach system. Syntax:
%
%                 summary_basis_opts(spin_system)
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

function summary_basis_opts(spin_system)

% Check consistency
grumble(spin_system);

% Report the formalism
switch spin_system.bas.formalism
    case 'zeeman-wavef'
        report(spin_system,'Zeeman basis set using wavefunction formalism.');
    case 'zeeman-hilb'
        report(spin_system,'Zeeman basis set using Hilbert space matrix formalism.');
    case 'zeeman-liouv'
        report(spin_system,'Zeeman basis set using Liouville space matrix formalism.');
    case 'sphten-liouv'
        report(spin_system,'spherical tensor basis set using Liouville space matrix formalism.');
    otherwise
        error('unrecognized formalism - see the basis preparation section of the manual.');
end

% Report the approximation
if strcmp(spin_system.bas.formalism,'sphten-liouv')
    switch spin_system.bas.approximation
        case 'IK-0'
            report(spin_system,['IK-0 approximation - all correlations of all spins up to order ' int2str(spin_system.bas.level)]);
        case 'IK-1'
            report(spin_system,['IK-1 approximation - spin correlations up to order ' int2str(spin_system.bas.level)...
                                ' between directly coupled spins.']);
            report(spin_system,['IK-1 approximation - spin correlations up to order ' int2str(spin_system.bas.space_level)...
                                ' between all spins within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
        case 'IK-2'
            report(spin_system, 'IK-2 approximation - spin correlations involving all nearest neighbours of each spin on the coupling graph.');
            report(spin_system,['IK-2 approximation - spin correlations up to order ' int2str(spin_system.bas.space_level)...
                                ' between all spins within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
        case 'IK-DNP'
            report(spin_system,['IK-DNP approximation | Max inter-electron correlation level:   ' int2str(spin_system.bas.level(1))]);
            report(spin_system,['IK-DNP approximation | max electron-nuclear correlation level: ' int2str(spin_system.bas.level(2))]);
            report(spin_system,['IK-DNP approximation | max inter-nuclear correlation level:    ' int2str(spin_system.bas.level(3))]);
            report(spin_system, 'IK-DNP approximation | with nearest neighbours on the coupling graph.');
        case 'none'
            report(spin_system, 'starting with complete basis set on all spins...');
        otherwise
            error('unrecognized approximation level - see the basis set preparation manual.');
    end
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

% According to a trade legend, Uhlenbeck and Goudsmit (students of
% Ehrenfest when they stumbled upon the concept of spin) presented
% it to Ehrenfest and said, in effect "here's our theory, but don't
% publish it - it can't be right". He submitted it anyway with the
% justification that they were "young enough to be able to afford
% a stupidity".

