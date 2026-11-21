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
%
% <https://spindynamics.org/wiki/index.php?title=summary_basis.m>

function summary_basis(spin_system)

% Check consistency
grumble(spin_system);

% Report every substance in its local spin columns
for s=1:spin_system.bas.nsubst

    % Identify the substance and its local dimensions
    spins=spin_system.chem.parts{s}; nspins=numel(spins);
    nstates=spin_system.bas.nstates(s);
    report(spin_system,['chemical substance ' num2str(s)]);

    % Print spherical tensor labels unless the table is too large
    if nstates>spin_system.tols.basis_hush
        report(spin_system,['over ' num2str(spin_system.tols.basis_hush) ...
                            ' states in the basis - printing suppressed.']);
    elseif strcmp(spin_system.bas.formalism,'sphten-liouv')
        report(spin_system,'final basis set summary (L,M quantum numbers in irreducible spherical tensor products).');

        % Size columns from local ranks, signed projections, and global indices
        [ranks,projs]=lin2lm(spin_system.bas.basis{s});
        rank_width=numel(num2str(max(ranks,[],'all')));
        proj_width=max(numel(num2str(min(projs,[],'all'))),...
                       numel(num2str(max(projs,[],'all'))));
        col_width=max([8 rank_width+proj_width+5 ...
                       numel(num2str(max(spins)))+2]);
        index_width=max(7,numel(num2str(spin_system.bas.offsets(s)+nstates))+1);

        % Retain the header inset used by the ordinary eight-character columns
        header_line=blanks(index_width+col_width*nspins); header_line(1)='N';
        for k=1:nspins
            spin_number=num2str(spins(k));
            header_line(index_width+1+col_width*(k-1)+(1:length(spin_number)))=spin_number;
        end
        report(spin_system,deblank(header_line));
        for n=1:nstates

            % Format the global state index and local tensor labels
            current_line=blanks(index_width+col_width*nspins);
            state_number=num2str(spin_system.bas.offsets(s)+n);
            current_line(1:length(state_number))=state_number;
            for k=1:nspins
                state_token=['(' num2str(ranks(n,k)) ',' num2str(projs(n,k)) ')'];
                current_line(index_width+col_width*(k-1)+(1:length(state_token)))=state_token;
            end
            report(spin_system,current_line);
        end
        report(spin_system,' ');
    end

    % Report the fraction of the full local state space
    full_dim=prod(spin_system.comp.mults(spins));
    if ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
        full_dim=full_dim^2;
    end
    report(spin_system,['state space dimension ' num2str(nstates) ...
                        ' (' num2str(100*nstates/full_dim) ...
                        '% of the full state space).']);

end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if isfield(spin_system.bas,'basis')&&~iscell(spin_system.bas.basis)
    error('Spinach:basis:retiredGlobalBasis',...
          'the global bas.basis matrix is retired; use bas.basis{n} and bas.offsets from basis().');
end
if isfield(spin_system.bas,'irrep')
    error('Spinach:basis:retiredIrrep',...
          'bas.irrep is retired; use bas.sym_fact(n).irr_projectors and irr_dimensions.');
end
end

% Linux is only free if your time has no value.
%
% Jamie Zawinski

