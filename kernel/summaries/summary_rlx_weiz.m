% Prints Weizmann DNP relaxation-rate summary for a Spinach system. Syntax:
%
%                 summary_rlx_weiz(spin_system)
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

function summary_rlx_weiz(spin_system)

% Check consistency
grumble(spin_system);

% Print the relaxation-rate table
report(spin_system,' ');
report(spin_system,'==== Weizmann DNP relaxation theory ====');
report(spin_system,['Electron R1: ' num2str(spin_system.rlx.weiz_r1e) ' Hz']);
report(spin_system,['Electron R2: ' num2str(spin_system.rlx.weiz_r2e) ' Hz']);
report(spin_system,['Nuclear R1:  ' num2str(spin_system.rlx.weiz_r1n) ' Hz']);
report(spin_system,['Nuclear R2:  ' num2str(spin_system.rlx.weiz_r2n) ' Hz']);
[rows,cols,vals]=find(spin_system.rlx.weiz_r1d);
for n=1:numel(vals)
    report(spin_system,['Inter-nuclear dipolar R1(' num2str(rows(n)) ',' num2str(cols(n)) '): ' num2str(vals(n)) ' Hz']);
end
[rows,cols,vals]=find(spin_system.rlx.weiz_r2d);
for n=1:numel(vals)
    report(spin_system,['Inter-nuclear dipolar R2(' num2str(rows(n)) ',' num2str(cols(n)) '): ' num2str(vals(n)) ' Hz']);
end
report(spin_system,'========================================');

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

