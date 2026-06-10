% Prints Nottingham DNP relaxation-rate summary for a Spinach system. Syntax:
%
%                 summary_rlx_nott(spin_system)
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

function summary_rlx_nott(spin_system)

% Check consistency
grumble(spin_system);

% Print the relaxation-rate table
report(spin_system,' ');
report(spin_system,'=== Nottingham DNP relaxation theory ===');
report(spin_system,['Electron R1: ' num2str(spin_system.rlx.nott_r1e) ' Hz']);
report(spin_system,['Electron R2: ' num2str(spin_system.rlx.nott_r2e) ' Hz']);
report(spin_system,['Nuclear R1:  ' num2str(spin_system.rlx.nott_r1n) ' Hz']);
report(spin_system,['Nuclear R2:  ' num2str(spin_system.rlx.nott_r2n) ' Hz']);
report(spin_system,'========================================');

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

% As a man, sometimes you have to make a choice
% between mocking astrology and getting laid.
%
% Internet wisdom

