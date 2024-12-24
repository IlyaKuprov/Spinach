% Emulates a strong homospoil pulse - only zero-frequency states with
% respect to the carrier frequencies (chemical shifts are not conside-
% red) survive the process. Syntax:
%
%               rho=homospoil(spin_system,rho,zqc_flag)
%
% Parameters:
%
%         rho - a state vector or a horizontal stack thereof
%
%    zqc_flag - a flag controlling the fate of zero-quantum 
%               coherences. If set to 'keep', causes ZQCs to
%               survive the process, approximating experimen-
%               tal behaviour. If set to 'destroy', wipes the
%               zero-quantum coherences - only the longitudi-
%               nal states survive the process. 
%
%               The flag is ignored in zeeman-hilb and zeeman-
%               liouv formalisms, where the effect is always
%               to destroy everything except the diagonal of
%               the density matrix.
%
% Outputs:
%
%        rho  - the state vector(s) with only the longitudi-
%               nal or only the zero-quantum states kept
%
% Note: this function is only available for sphten-liouv formalism; it
%       supports Fokker-Planck direct products.
%
% Note: this is a purely mathematical filter that only mimics - in an
%       idealised way - the effect of a real homospoil pulse. Essenti-
%       ally, it searches the density matrix for any transverse state
%       populations and zeroes them out. If the flag is set, zero-qua-
%       ntum coherences are also erased.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=homospoil.m>

function rho=homospoil(spin_system,rho,zqc_flag)

% Check consistency
grumble(rho,zqc_flag);

% In Hilbert space, only keep the diagonal
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    dim=size(rho,1); rhod=diag(rho);
    rho=spdiags(rhod,0,dim,dim); return;
end

% In Zeeman basis of Liouville space, fold back
% into Hilbert space, keep diagonal, and unfold
if strcmp(spin_system.bas.formalism,'zeeman-liouv')
    dim=sqrt(size(rho,1));
    rho=reshape(rho,[dim dim]);
    rho=spdiags(diag(rho),0,dim,dim);
    rho=rho(:); return;
end

% Store dimension statistics
spn_dim=size(spin_system.bas.basis,1);
spc_dim=numel(rho)/spn_dim;
problem_dims=size(rho);

% Fold indirect dimensions
rho=reshape(rho,[spn_dim spc_dim]);

% Pull the projection information from the basis
[~,M]=lin2lm(spin_system.bas.basis);

% Filter the state vector
switch zqc_flag
    
    case 'keep'
        
        % Find the states that have zero carrier frequency and kill everything else
        rho(abs(sum(repmat(spin_system.inter.basefrqs,size(spin_system.bas.basis,1),1).*M,2))>1e-6,:)=0;
    
    case 'destroy'
        
        % Find the longitudinal states and kill everything else
        rho(sum(abs(M),2)>0,:)=0;
    
    otherwise
        
        % Complain and bomb out
        error('unknown ZQC flag.');
        
end

% Unfold indirect dimensions
rho=reshape(rho,problem_dims);

% Report overly destructive calls
if norm(rho,1)<1e-10
    report(spin_system,'WARNING - all magnetization appears to have been destroyed by this call.');
end

end

% Consistency enforcement
function grumble(rho,zqc_flag)
if ~isnumeric(rho)
    error('the state vector(s) must be numeric.');
end
if ~ischar(zqc_flag)
    error('zqc_flag parameter must be a character string.');
end
if ~ismember(zqc_flag,{'keep','destroy'})
    error('the available values for zqc_flag are ''keep'' and ''destroy''.');
end
end

% Rocket science has been mythologized all out of proportion to
% its true difficulty.
%
% John Carmack

