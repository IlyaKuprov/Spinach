% Point-to-point magnetisation flux reporting. Syntax:
%
%                    summary_rlx_flux(spin_system)
%
% Parameters:
%
%    spin_system - Spinach data structure with validated flux inputs
%
% Outputs:
%
%    This function prints flux rates to the configured output stream.
%
% i.kuprov@soton.ac.uk

function summary_rlx_flux(spin_system)

% Check consistency
grumble(spin_system);

% Report flux inputs before basis construction
if isfield(spin_system.chem,'flux_rate')&&nnz(spin_system.chem.flux_rate)>0
    switch spin_system.chem.flux_type
        case 'intramolecular'
            report(spin_system,'intramolecular magnetisation flux, correlations will be kept.');
        case 'intermolecular'
            report(spin_system,'intermolecular magnetisation flux, correlations will be damped.');
    end
    [rows,cols,vals]=find(spin_system.chem.flux_rate);
    report(spin_system,'point-to-point flux rates:');
    report(spin_system,'===============================');
    report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
    report(spin_system,'-------------------------------');
    for n=1:numel(vals)
        report(spin_system,sprintf(' %-3d      %-3d      %+0.3e',rows(n),cols(n),vals(n)));
    end
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

