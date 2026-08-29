% Obliterates all interactions and populations in the subspace of states
% that involve the specified spins in any way. The specified spins would
% not contribute to the system dynamics until the Liouvillian is rebuilt
% from scratch. Syntax:
%
%                 [L,rho]=decouple(spin_system,L,rho,spins)
%
% Parameters:
%
%     L     - Liouvillian superoperator or, in zeeman-hilb,
%             the Hamiltonian; this may be left empty
%
%     rho   - state vector or a horizontal stack thereof or,
%             in zeeman-hilb, a density matrix or a horizon-
%             tal stack thereof; this may be left empty
%
%     spins - spins to be wiped, specified either by name, e.g.
%             {'13C','1H'}, or by a list of numbers, e.g. [1 2]
%
% Outputs:
%
%     rho   - state vector(s) with all populations of the 
%             states involving the target spins set to zero
%
%     L     - Liouvillian superoperator with all terms in-
%             volving the target spins set to zero 
%
% Note: this function is an analytical equivalent of a perfect decoup-
%       ling pulse sequence on the specified spins.
%
% Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
%       hilb formalism; Fokker-Planck direct products are supported in
%       the Liouville space formalisms. In the Zeeman formalisms the
%       operation is an exact projection onto the subspace where the
%       decoupled spins carry only their identity component, because
%       spin involvement is not diagonal in the Zeeman basis. In
%       zeeman-hilb, the Hamiltonian and the density matrices are
%       stretched into Liouville space, projected there, and folded
%       back; this replaces every decoupled-spin factor of the Hamil-
%       tonian by its identity component average.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=decouple.m>

function [L,rho]=decouple(spin_system,L,rho,spins)

% Return if the spin list is empty
if isempty(spins), return; end

% Check consistency
grumble(spin_system,L,rho,spins);

% Find the spins to be decoupled
if isnumeric(spins)

    % Find spins by numbers
    dec_mask=false(1,spin_system.comp.nspins); dec_mask(spins)=true;

else

    % Find spins by name
    dec_mask=ismember(spin_system.comp.isotopes,spins);

end

% Inform the user
report(spin_system,[num2str(nnz(dec_mask)) ' spins to be frozen and depopulated.']);

% Get the spin space dimension
spn_dim=size(spin_system.bas.basis,1);

% Build the wipeout machinery
switch spin_system.bas.formalism

    case 'sphten-liouv'

        % Get the list of states to be wiped
        zero_mask=(sum(spin_system.bas.basis(:,dec_mask),2)~=0);

    case {'zeeman-liouv','zeeman-hilb'}

        % Hilbert space dimension and multiplicities
        dim=prod(spin_system.comp.mults);
        mults=spin_system.comp.mults;

        % Identity component projector of the decoupled spins
        P=speye(dim^2);
        for n=find(dec_mask)
            idch=sparse(dim^2,dim^2);
            for p=1:mults(n)
                for q=1:mults(n)
                    A=kron(kron(speye(prod(mults(1:(n-1)))),...
                                sparse(p,q,1,mults(n),mults(n))),...
                           speye(prod(mults((n+1):end))));
                    idch=idch+kron(conj(A),A)/mults(n);
                end
            end
            P=P*idch;
        end

end

% Process the Liouvillian
if (nargout>0)&&(~isempty(L))

    % Get dimension statistics
    spc_dim=size(L,1)/spn_dim;

    % Inform the user
    report(spin_system,['space sub-problem dimension: ' num2str(spc_dim)]);
    report(spin_system,['spin sub-problem dimension:  ' num2str(spn_dim)]);

    % Apply the wipeout machinery
    switch spin_system.bas.formalism

        case 'sphten-liouv'

            % Kron the list into the Fokker-Planck space
            fp_zero_mask=logical(kron(ones(spc_dim,1),zero_mask));

            % Inform the user
            report(spin_system,['zeroing ' num2str(nnz(fp_zero_mask))...
                                ' rows and columns in the Liouvillian.']);

            % Apply the zero mask
            L(fp_zero_mask,:)=0; L(:,fp_zero_mask)=0;

        case 'zeeman-liouv'

            % Kron the projector into the Fokker-Planck space
            fp_proj=kron(speye(spc_dim),P);

            % Apply the projector from both sides
            L=fp_proj*L*fp_proj;

        case 'zeeman-hilb'

            % Stretch, project, and fold the Hamiltonian
            L=reshape(P*L(:),size(L));

    end

    % Re-evaluate sparsity
    L=clean_up(spin_system,L,spin_system.tols.liouv_zero);

end

% Process state vector stack
if (nargout>1)&&(~isempty(rho))

    % Get dimension statistics
    spc_dim=size(rho,1)/spn_dim;

    % Inform the user
    report(spin_system,['space sub-problem dimension: ' num2str(spc_dim)]);
    report(spin_system,['spin sub-problem dimension:  ' num2str(spn_dim)]);

    % Apply the wipeout machinery
    switch spin_system.bas.formalism

        case 'sphten-liouv'

            % Kron the list into the Fokker-Planck space
            fp_zero_mask=logical(kron(ones(spc_dim,1),zero_mask));

            % Inform the user
            report(spin_system,['zeroing ' num2str(nnz(fp_zero_mask))...
                                ' rows in the state vector.']);

            % Apply the zero mask
            rho(fp_zero_mask,:)=0;

        case 'zeeman-liouv'

            % Kron the projector into the Fokker-Planck space
            fp_proj=kron(speye(spc_dim),P);

            % Apply the projector
            rho=fp_proj*rho;

        case 'zeeman-hilb'

            % Stretch, project, and fold the state stack
            rho=reshape(P*reshape(rho,spn_dim^2,[]),size(rho));

    end

end

end

% Consistency enforcement
function grumble(spin_system,L,rho,spins)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})
    error('analytical decoupling is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.');
end
if (~isempty(rho))&&(~isempty(L))&&(size(L,2)~=size(rho,1))
    error('matrix dimensions of L and rho must agree.');
end
if size(L,1)~=size(L,2)
    error('Liouvillian must be a square matrix.');
end
if (~isnumeric(spins))&&(~iscell(spins))
    error('spins parameter must either be a list of numbers or a cell array of strings.');
end
if iscell(spins)&&any(~ismember(spins,spin_system.comp.isotopes))
    error('the system does not contain the spins specified.');
end
if isnumeric(spins)&&(size(spins,1)~=1)
    error('if spins are specified by number, a row vector of numbers must be used.');
end
if isnumeric(spins)&&(max(spins)>spin_system.comp.nspins)
    error('the spin number specified is greater than the number of spins in the system.');
end
if isnumeric(spins)&&(any(~isreal(spins))||any(spins<1)||any(mod(spins,1)~=0))
    error('spin numbers must be real positive integers.');
end
end

% It's not worth doing something unless you were doing something that
% someone, somewere, would much rather you weren't doing.
%
% Terry Pratchett 

