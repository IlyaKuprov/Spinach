% Moves a zeeman-hilb simulation context into Liouville space. When
% the formalism specified in the spin system object is 'zeeman-hilb',
% this function projects the evolution generators into Liouville
% space, converts the standard state-like and operator-like fields
% of the parameters structure, rebuilds the block dimensions, mig-
% rates the symmetry irrep projectors into the adjoint representa-
% tion, and sets the formalism to 'zeeman-liouv'; for all other
% formalisms, every argument is returned unchanged. This makes
% Liouville-space pulse sequences callable with zeeman-hilb
% inputs. Syntax:
%
%          [spin_system,parameters,H,R,K]=...
%          sim2liouv(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    spin_system  - Spinach spin system object
%
%    parameters   - pulse sequence parameters structure; the
%                   state-like fields rho0, coil, and screen
%                   (matrices or their horizontal concatena-
%                   tions) are stretched into state vectors,
%                   and the operator-like fields pulse_op,
%                   mw_oper, ez_oper, and homodec_oper become
%                   commutation superoperators, when present
%
%    H            - Hamiltonian operator, converted into a
%                   commutation superoperator; an empty
%                   matrix is passed through
%
%    R            - relaxation matrix, converted into an
%                   anticommutation superoperator with the
%                   unit state exempted from damping; an
%                   empty matrix is passed through
%
%    K            - kinetics matrix, converted into an
%                   anticommutation superoperator; an empty
%                   matrix is passed through
%
% Outputs:
%
%    spin_system  - spin system object with zeeman-liouv
%                   formalism and basis information
%
%    parameters   - parameters structure with the standard
%                   fields converted into Liouville space
%
%    H,R,K        - Liouville space evolution generators
%
% Note: a Hilbert space density matrix block S(n)*Y*S(k)' maps to
%       the state vector kron(conj(S(k)),S(n))*Y(:), and so each
%       ordered pair of Hilbert space irrep projectors yields the
%       Liouville space irrep projector kron(conj(S(k)),S(n)).
%       Every such subspace is invariant under superoperators
%       built from symmetry-respecting Hilbert space generators;
%       unpopulated subspaces are dropped by reduce.m at run time
%       in the usual way.
%
% Note: the anticommutation superoperator of a Hilbert space rela-
%       xation matrix damps the unit state, which the Liouville
%       space branch of relaxation.m never does. The row and the
%       column of the unit state are therefore projected out of
%       the converted R, so that the unit state is neither damped
%       nor a source of relaxation and the trace is conserved. For
%       the scalar damping matrix that the kernel builds in Hilbert
%       space the result coincides with the Liouville space damp
%       operator. The unit state spans all diagonal irrep pairs,
%       and those are merged into one subspace so that the projec-
%       ted R stays block-diagonal in the irrep table that reduce.m
%       evolves independently.
%
% Note: segmented Hilbert inputs must be block diagonal in substance;
%       cross-substance entries in generators and parameter matrices
%       are rejected rather than discarded during conversion.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sim2liouv.m>

function [spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Only the zeeman-hilb formalism needs the move
if strcmp(spin_system.bas.formalism,'zeeman-hilb')

    % Inform the user
    report(spin_system,'projecting zeeman-hilb simulation into Liouville space...');

    % Convert generators using explicit Hilbert space substance blocks
    dims=spin_system.bas.nstates; hdim=sum(dims);
    generators={H,R,K}; types={'comm','acomm','acomm'};
    for n=1:numel(generators)
        if ~isempty(generators{n})
            blocks=mat2cell(generators{n},dims,dims);
            generators{n}=hilb2liouv(blocks(1:spin_system.bas.nsubst+1:end),types{n});
        end
    end
    [H,R,K]=generators{:};

    % Stretch each state-like parameter within each substance
    fields={'rho0','coil','screen'};
    for n=1:numel(fields)
        if isfield(parameters,fields{n})
            states=parameters.(fields{n});
            blocks=cell(spin_system.bas.nsubst,1);
            for k=1:spin_system.bas.nsubst
                idx=spin_system.bas.offsets(k)+(1:dims(k));
                cols=reshape(idx(:)+hdim*(0:size(states,2)/hdim-1),1,[]);
                blocks{k}=reshape(states(idx,cols),dims(k)^2,[]);
            end
            parameters.(fields{n})=vertcat(blocks{:});
        end
    end

    % Convert operator-like parameters within each substance
    fields={'pulse_op','mw_oper','ez_oper','homodec_oper'};
    for n=1:numel(fields)
        if isfield(parameters,fields{n})
            blocks=mat2cell(parameters.(fields{n}),dims,dims);
            parameters.(fields{n})=hilb2liouv(blocks(1:spin_system.bas.nsubst+1:end),'comm');
        end
    end

    % Compile the Liouville dimensions and refresh the cache identity
    spin_system.bas.nstates=dims.^2;
    spin_system.bas.offsets=[0;cumsum(spin_system.bas.nstates)];
    spin_system.bas.basis_hash=md5_hash({spin_system.bas.basis,...
                                       spin_system.bas.nstates,spin_system.chem.parts});

    % Migrate the irreps into the adjoint representation
    for s=1:spin_system.bas.nsubst

        % Grab the Hilbert space irreps
        hs_irreps=spin_system.bas.sym_fact(s);
        n_irreps=numel(hs_irreps.irr_dimensions);

        % Preallocate the Liouville space irrep array
        ls_irreps=repmat(struct('projector',[],'dimension',[]),n_irreps^2-n_irreps+1,1);

        % Diagonal irrep pairs share the unit state and are merged into the first subspace
        ls_irreps(1).projector=[]; ls_irreps(1).dimension=0; pair_idx=1;

        % Loop over ordered pairs of Hilbert space irreps
        for n=1:n_irreps
            for k=1:n_irreps

                % Build the irrep pair projector and dimension
                pair_proj=kron(conj(hs_irreps.irr_projectors{n}),hs_irreps.irr_projectors{k});
                pair_dim=hs_irreps.irr_dimensions(n)*hs_irreps.irr_dimensions(k);

                % Merge diagonal pairs, store off-diagonal pairs separately
                if n==k
                    ls_irreps(1).projector=[ls_irreps(1).projector pair_proj];
                    ls_irreps(1).dimension=ls_irreps(1).dimension+pair_dim;
                else
                    pair_idx=pair_idx+1;
                    ls_irreps(pair_idx).projector=pair_proj;
                    ls_irreps(pair_idx).dimension=pair_dim;
                end

            end
        end

        % Write the Liouville space irreps
        spin_system.bas.sym_fact(s).irr_dimensions=[ls_irreps.dimension]';
        spin_system.bas.sym_fact(s).irr_projectors={ls_irreps.projector};
        report(spin_system,['Hilbert space irreps migrated into Liouville space, '...
                            num2str(numel(ls_irreps)) ' irrep pair subspaces.']);

    end

    % Update the formalism setting
    spin_system.bas.formalism='zeeman-liouv';

    % Project the unit state out of the relaxation superoperator
    if ~isempty(R)
        units=cell(spin_system.bas.nsubst,1);
        for n=1:spin_system.bas.nsubst
            unit=speye(dims(n)); units{n}=unit(:)/sqrt(dims(n));
        end
        U=blkdiag(units{:});
        R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U';
        report(spin_system,'unit state exempted from the projected relaxation superoperator.');
    end

end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv',...
                                        'zeeman-hilb','zeeman-wavef'})
    error('unrecognised formalism in spin_system.bas.formalism.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))
    error('H, R, and K must be numeric arrays.');
end
if strcmp(spin_system.bas.formalism,'zeeman-hilb')&&spin_system.bas.nsubst>1
    fields={'pulse_op','mw_oper','ez_oper','homodec_oper','rho0','coil','screen'};
    fields=fields(isfield(parameters,fields));
    matrices=[{H,R,K} cell(1,numel(fields))];
    for n=1:numel(fields)
        matrices{n+3}=parameters.(fields{n});
    end
    membership=repelem((1:spin_system.bas.nsubst)',spin_system.bas.nstates);
    for n=1:numel(matrices)
        [rows,cols]=find(matrices{n});
        cols=mod(cols-1,spin_system.bas.offsets(end))+1;
        if any(membership(rows)~=membership(cols))
            error('Spinach:sim2liouv:crossSubstance',...
                  'Hilbert inputs must not contain cross-substance matrix entries.');
        end
    end
end
end

% To achieve great things, two things are needed: a plan,
% and not quite enough time.
%
% Leonard Bernstein

