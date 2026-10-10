% Returns a matrix that converts state vectors written in the 
% spherical tensor basis set used by Spinach into state vectors
% written in the Zeeman basis set in Liouville space. Syntax:
%
%                  P=sphten2zeeman(spin_system)
%
% Parameters:
%
%    spin_system - main Spinach data structure using
%                  sphten-liouv formalism and inclu-
%                  ding basis set information
%
% Outputs:
%
%    P - projector matrix that is to be used in the fol-
%        lowing way:
%
%                   rho_zeeman=P*rho_sphten
%
% Note: the projector need not be square and may be huge. Each
%       substance is converted separately, with destination dimension
%       D_n^2. The unit coordinate maps to vec(I_D_n), so the source
%       unit coordinate equals the Hilbert trace divided by D_n.
%       For trace-equals-concentration CWDM coordinates, divide each
%       destination substance block of P by its local D_n explicitly.
%
% ilya.kuprov@weizmann.ac.il
% enu.jamila@proton.me
%
% <https://spindynamics.org/wiki/index.php?title=sphten2zeeman.m>

function P=sphten2zeeman(spin_system)

% Check consistency
grumble(spin_system);

% Build one conversion matrix per substance
blocks=cell(spin_system.bas.nsubst,1);
for s=1:spin_system.bas.nsubst

    % Read the local descriptor and multiplicities
    descriptor=spin_system.bas.basis{s};
    mults=spin_system.comp.mults(spin_system.chem.parts{s});
    destin_norm=sqrt(prod(mults));
    block=spalloc(prod(mults.^2),spin_system.bas.nstates(s),0);

    % Convert the local spherical tensor basis
    parfor n=1:size(descriptor,1)

        % Form the tensor product for this state
        rho=sparse(1);
        for k=1:size(descriptor,2)
            ists=irr_sph_ten(mults(k)); %#ok<PFBNS>
            rho=kron(rho,ists{descriptor(n,k)+1});
        end

        % Preserve the source and destination normalisations
        source_norm=norm(rho(:),2);
        block(:,n)=destin_norm*rho(:)/source_norm; %#ok<SPRIX>

    end

    % Store the local conversion
    blocks{s}=block;

end

% Never introduce inter-substance coherences
P=blkdiag(blocks{:});

end

% Consistency enforcement
function grumble(spin_system)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function is only available for sphten-liouv formalism.');
end
end

% Nurture your minds with great thoughts. To believe
% in the heroic makes heroes.
%
% Benjamin Disraeli

