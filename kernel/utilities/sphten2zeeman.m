% Returns a matrix that converts state vectors written in the 
% spherical tensor basis set used by Spinach into state vectors
% written in the Zeeman basis set in Liouville space. Syntax:
%
%                  P=sphten2zeeman(spin_system)
%
% Outputs:
%
%    P - projector matrix that is to be used in the fol-
%        lowing way:
%
%                   rho_zeeman=P*rho_sphten
%
% Note: the matrix need not be square and may be huge.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sphten2zeeman.m>

function P=sphten2zeeman(spin_system)

% Check consistency
grumble(spin_system);

% Preallocate the answer
P=spalloc(prod(spin_system.comp.mults.^2),size(spin_system.bas.basis,1),0);

% Loop over the basis set
for n=1:size(spin_system.bas.basis,1)

    % Get the state going
    rho=sparse(1);
    
    % Loop over the elements
    for k=1:size(spin_system.bas.basis,2)
        
        % Get the spherical tensors for the current spin
        ists=irr_sph_ten(spin_system.comp.mults(k));
        
        % Multiply into the state
        rho=kron(rho,ists{spin_system.bas.basis(n,k)+1});
        
    end
    
    % Stretch and write down
    P(:,n)=rho(:)/norm(rho(:),2); %#ok<SPRIX>
    
end

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

