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
                        rate_text(spin_system.rlx.r1_rates{n}) '  '...
                        rate_text(spin_system.rlx.r2_rates{n}) '  '...
                        spin_system.comp.labels{n}]);
end
report(spin_system,'========================================');

end

% Formats a scalar, tensor, or function-handle relaxation rate
function rate_str=rate_text(rate)
if isnumeric(rate)&&isscalar(rate)
    rate_str=num2str(rate,'%+0.5e   ');
elseif isnumeric(rate)
    rate_str='anisotropic    ';
else
    rate_str='orientation    ';
end
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

