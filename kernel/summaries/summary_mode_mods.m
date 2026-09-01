% Prints the summary of spin Hamiltonian modulation by bosonic
% mode coordinates: derivatives of spin-spin coupling tensors and
% of effective local fields with respect to dimensionless mode
% displacement coordinates. Syntax:
%
%               summary_mode_mods(spin_system,header)
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
% <https://spindynamics.org/wiki/index.php?title=summary_mode_mods.m>

function summary_mode_mods(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the summary table
report(spin_system,header);
report(spin_system,'=============================================================');
report(spin_system,'Mode A  Mode B  Order  Modulation target    Norm, Hz        ');
report(spin_system,'-------------------------------------------------------------');

% Print coupling tensor derivatives
[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.coupling_mod));
for n=1:numel(rows)
    deriv_orders=spin_system.inter.modes.coupling_mod{rows(n),cols(n)};
    for m=1:numel(deriv_orders)
        if ~isempty(deriv_orders{m})
            [spr,spc]=find(~cellfun(@isempty,deriv_orders{m}));
            for k=1:numel(spr)
                target=['coupling {' num2str(spr(k)) ',' num2str(spc(k)) '}'];
                report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8)...
                                    pad(num2str(m),7) pad(target,21)...
                                    num2str(norm(deriv_orders{m}{spr(k),spc(k)},'fro')/(2*pi),'%+10.4e')]);
            end
        end
    end
end

% Print effective field derivatives
[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.zeeman_mod));
for n=1:numel(rows)
    deriv_orders=spin_system.inter.modes.zeeman_mod{rows(n),cols(n)};
    for m=1:numel(deriv_orders)
        if ~isempty(deriv_orders{m})
            spr=find(~cellfun(@isempty,deriv_orders{m}));
            for k=1:numel(spr)
                target=['field, spin ' num2str(spr(k))];
                report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8)...
                                    pad(num2str(m),7) pad(target,21)...
                                    num2str(norm(deriv_orders{m}{spr(k)},2)/(2*pi),'%+10.4e')]);
            end
        end
    end
end
report(spin_system,'=============================================================');

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

