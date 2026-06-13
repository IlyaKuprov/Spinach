% Prints T1 and T2 relaxation-rate summary for a Spinach system. Syntax:
%
%                 summary_rlx_t1_t2(spin_system,header)
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
% <https://spindynamics.org/wiki/index.php?title=summary_rlx_t1_t2.m>

function summary_rlx_t1_t2(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the relaxation-rate table
report(spin_system,header);
report(spin_system,'========================================');
report(spin_system,'N    Spin        R1             R2      ');
report(spin_system,'----------------------------------------');
for n=1:spin_system.comp.nspins
    report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                        strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                        num2str(spin_system.rlx.r1_rates(n),'%+0.5e   ') '  '...
                        num2str(spin_system.rlx.r2_rates(n),'%+0.5e   ') '  '...
                        spin_system.comp.labels{n}]);
end
report(spin_system,'========================================');

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

% Whether you are more afraid of the forces of 
% order or the forces of chaos is generally a
% matter of disposition.
%
% Lionel Shriver

