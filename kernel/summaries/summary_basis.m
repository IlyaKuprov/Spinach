% Prints basis-set state summary for a Spinach system. Syntax:
%
%                 summary_basis(spin_system)
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

function summary_basis(spin_system)

% Check consistency
grumble(spin_system);

% Get the basis dimension
nstates=size(spin_system.bas.basis,1);
if nstates > spin_system.tols.basis_hush
    report(spin_system,['over ' num2str(spin_system.tols.basis_hush) ' states in the basis - printing suppressed.']);
else
    report(spin_system,'final basis set summary (L,M quantum numbers in irreducible spherical tensor products).')
    report(spin_system,['N       ' blanks(length(num2str(spin_system.comp.nspins))) ...
                        num2str(1:spin_system.comp.nspins,['%d ' blanks(7-length(num2str(spin_system.comp.nspins)))])]);
    for n=1:nstates
        current_line=blanks(7+8*spin_system.comp.nspins); spin_number=num2str(n);
        current_line(1:length(spin_number))=spin_number;
        for k=1:spin_system.comp.nspins
            [L,M]=lin2lm(spin_system.bas.basis(n,k));
            current_line(7+8*(k-1)+1)='(';
            current_line(7+8*(k-1)+2)=num2str(L);
            current_line(7+8*(k-1)+3)=',';
            proj=num2str(M);
            switch length(proj)
                case 1
                    current_line(7+8*(k-1)+4)=proj;
                    current_line(7+8*(k-1)+5)=')';
                case 2
                    current_line(7+8*(k-1)+4)=proj(1);
                    current_line(7+8*(k-1)+5)=proj(2);
                    current_line(7+8*(k-1)+6)=')';
            end
        end
        report(spin_system,current_line);
    end
    report(spin_system,' ');
end
report(spin_system,['state space dimension ' num2str(nstates) ...
                    ' (' num2str(100*nstates/(prod(spin_system.comp.mults)^2))...
                    '% of the full state space).']);

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

