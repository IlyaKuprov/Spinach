% Moves a zeeman-hilb simulation context into Liouville space. When
% the formalism specified in the spin system object is 'zeeman-hilb',
% this function projects the evolution generators into Liouville
% space, converts the standard state-like and operator-like fields
% of the parameters structure, rebuilds the basis index table, mig-
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
%                   mw_oper, and ez_oper are converted into
%                   commutation superoperators, when present
%
%    H            - Hamiltonian operator, converted into a
%                   commutation superoperator; an empty
%                   matrix is passed through
%
%    R            - relaxation matrix, converted into an
%                   anticommutation superoperator; an empty
%                   matrix is passed through
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

    % Project the evolution generators into Liouville space
    H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm');

    % Get the Hilbert space basis table and dimension
    zbas=spin_system.bas.basis; hdim=size(zbas,1);

    % Stretch the state-like parameters
    if isfield(parameters,'rho0')
        parameters.rho0=reshape(parameters.rho0,hdim^2,[]);
    end
    if isfield(parameters,'coil')
        parameters.coil=reshape(parameters.coil,hdim^2,[]);
    end
    if isfield(parameters,'screen')
        parameters.screen=reshape(parameters.screen,hdim^2,[]);
    end

    % Project the operator-like parameters
    if isfield(parameters,'pulse_op')
        parameters.pulse_op=hilb2liouv(parameters.pulse_op,'comm');
    end
    if isfield(parameters,'mw_oper')
        parameters.mw_oper=hilb2liouv(parameters.mw_oper,'comm');
    end
    if isfield(parameters,'ez_oper')
        parameters.ez_oper=hilb2liouv(parameters.ez_oper,'comm');
    end

    % Rebuild the basis index table for the Liouville space
    spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))];

    % Migrate the irreps into the adjoint representation
    if isfield(spin_system.bas,'irrep')

        % Grab the Hilbert space irreps
        hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps);

        % Preallocate the Liouville space irrep array
        ls_irreps(n_irreps^2)=struct('projector',[],'dimension',[]);

        % Loop over ordered pairs of Hilbert space irreps
        for n=1:n_irreps
            for k=1:n_irreps

                % Build the irrep pair projector and dimension
                pair_idx=n_irreps*(n-1)+k;
                ls_irreps(pair_idx).projector=kron(conj(hs_irreps(n).projector),...
                                                        hs_irreps(k).projector);
                ls_irreps(pair_idx).dimension=hs_irreps(n).dimension*...
                                              hs_irreps(k).dimension;

            end
        end

        % Write the Liouville space irreps
        spin_system.bas.irrep=ls_irreps;
        report(spin_system,['Hilbert space irreps migrated into Liouville space, '...
                            num2str(n_irreps^2) ' irrep pair subspaces.']);

    end

    % Update the formalism setting
    spin_system.bas.formalism='zeeman-liouv';

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
end

% To achieve great things, two things are needed: a plan,
% and not quite enough time.
%
% Leonard Bernstein

