% Fast approximate simulation of isotropic continuous-wave ENDOR 
% pulse sequence - essentially an NMR spectrum weighted by hyper-
% fine couplings is recorded. Syntax:
%
%            fid=endor_cw(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   parameters.sweep         nuclear frequency sweep width, Hz
%
%   parameters.npoints       number of FID points to be computed
%
%   H  - Hamiltonian matrix, received from context function
%
%   R  - relaxation superoperator, received from context function
%
%   K  - kinetics superoperator, received from context function
%
% Outputs:
%
%   fid   - free induction decay whose Fourier transform 
%           approximates a CW ENDOR spectrum
%
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=endor_cw.m>

function fid=endor_cw(spin_system,parameters,H,R,K)

% Move into adjoint representation if needed
if strcmp(spin_system.bas.formalism,'zeeman-hilb')

    % Inform the user
    report(spin_system,'projecting zeeman-hilb simulation into Liouville space...');

    % Project into Liouville space
    H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm');

    % Stretch the state parameters
    if isfield(parameters,'rho0')
        parameters.rho0=reshape(parameters.rho0,size(spin_system.bas.basis,1)^2,[]);
    end
    if isfield(parameters,'coil')
        parameters.coil=reshape(parameters.coil,size(spin_system.bas.basis,1)^2,[]);
    end

    % Rebuild the basis index table for the Liouville space
    zbas=spin_system.bas.basis; hdim=size(zbas,1);
    spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))];

    % Discard Hilbert space symmetry data
    if isfield(spin_system.bas,'irrep')
        spin_system.bas=rmfield(spin_system.bas,'irrep');
        report(spin_system,'Hilbert space irreps discarded, proceeding without symmetry.');
    end

    % Update the formalism setting
    spin_system.bas.formalism='zeeman-liouv';

end

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;
        
% Nuclear pulse operator
Sp=operator(spin_system,'L+','nuclei'); Sy=(Sp-Sp')/2i;

% The initial state is a sum of nuclear Lz states weighted by HFCs
rho=zeros([size(L,1) 1],'like',1i);
for n=find(cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes))
    for k=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes))
        amplitude=trace(spin_system.inter.coupling.matrix{n,k})/3+...
                  trace(spin_system.inter.coupling.matrix{k,n})/3;
        rho=rho+abs(amplitude)*state(spin_system,{'Lz'},{k});
    end
end
rho=rho/norm(rho,2);

% Detect the nuclei
coil=state(spin_system,'L+','nuclei','cheap');
        
% Evolution time step
timestep=1/parameters.sweep;

% Pulse the nuclei
rho=step(spin_system,Sy,rho,pi/2);

% Acquire on the nuclei
fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');

% Frequency symmetrization
fid=real(fid);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('endor_cw: sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('endor_cw: parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('endor_cw: number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=1
    error('endor_cw: parameters.npoints array should have exactly one element.');
end
end

% It was the greatest sensation of existence: 
% not to trust, but to know.
%
% Ayn Rand, "Atlas Shrugged"

