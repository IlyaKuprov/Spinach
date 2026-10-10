% Slow passage detection - calculates spectrum values at the user-
% specified frequency positions using the Fourier transform of the
% Liouville - von Neumann equation. The biggest advantage over the
% fid+fft style detection is easy parallelization and the possibi-
% lity of getting spectrum values at specific frequencies without
% recalculating the entire free induction decay. Syntax:
%
%         spectrum=slowpass(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep         vector with two elements giving
%                             the spectrum frequency extents
%                             in Hz
%
%    parameters.npoints       number of points in the spectrum
%
%    parameters.rho0          initial state
%
%    parameters.coil          detection state
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    spectrum  - the spectrum of the system with the specified
%                starting state detected on the specified coil
%                state within the frequency interval requested
%
% Note: relaxation must be present in the system dynamics, or the 
%       matrix inversion operation would fail to converge. The re-
%       laxation matrix R must *not* be thermalized.
%
% Note: Liouville-space identity components are excluded only when
%       the identity sector is decoupled from spin order in both
%       directions, to within tols.liouv_zero. Its resolvent block
%       is shifted by 1 inverse second to remove stationary poles;
%       the spin-order block is unchanged. Spatial contexts supply
%       parameters.spc_dim for the space-times-spin embedding.
%       Coupled identity sectors, including selective reactions,
%       and wavefunction inputs retain their original resolvent.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=slowpass.m>

function spectrum=slowpass(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Get the frequency grid
freq_grid=2*pi*linspace(parameters.sweep(1),...
                        parameters.sweep(2),...
                        parameters.npoints)';

% Preallocate the answer
spectrum=zeros(size(freq_grid),'like',1i);

% Move into adjoint representation if needed
[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K);

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Compute subspace projectors
projectors=reduce(spin_system,L,parameters.coil);

% Get identity directions only in Liouville-space formalisms
U=sparse(size(L,1),0);
if ~strcmp(spin_system.bas.formalism,'zeeman-wavef')

    % Get normalised unit states, one per substance
    U=cell(1,spin_system.bas.nsubst);
    for n=1:spin_system.bas.nsubst
        unit_system=spin_system; unit_system.chem.concs(:)=0;
        unit_system.chem.concs(n)=1; U{n}=unit_state(unit_system);
    end
    U=[U{:}];

    % Embed the spin identities in every spatial basis coordinate
    if isfield(parameters,'spc_dim')
        U=kron(speye(parameters.spc_dim),U);
    end
end

% Loop over subspaces
for k=1:numel(projectors)
    
    % Project into the current subspace
    rho0_subs=projectors{k}'*parameters.rho0;
    coil_subs=projectors{k}'*parameters.coil;
    L_subs=projectors{k}'*L*projectors{k};
    Id_subs=speye(size(L_subs));

    % Normalise the projected identity directions
    U_subs=projectors{k}'*U; U_subs=U_subs(:,any(U_subs,1));
    U_subs=U_subs./sqrt(sum(abs(U_subs).^2,1));

    % Remove identities only when they decouple from spin order in both directions
    unit_block=U_subs'*L_subs*U_subs;
    if norm(L_subs*U_subs-U_subs*unit_block,1)<=spin_system.tols.liouv_zero&&...
       norm(U_subs'*L_subs-unit_block*U_subs',1)<=spin_system.tols.liouv_zero
        rho0_subs=rho0_subs-U_subs*(U_subs'*rho0_subs);
        coil_subs=coil_subs-U_subs*(U_subs'*coil_subs);
        L_subs=L_subs-1i*(U_subs*U_subs');
    end
    
    % Run backslash on the GPU if instructed
    if ismember('gpu',spin_system.sys.enable)
        
        % Inform the user
        report(spin_system,'using GPU backslash path...');
        
        % Move the objects to GPU
        rho0_subs=gpuArray(full(rho0_subs)); L_subs=gpuArray(full(L_subs));
        coil_subs=gpuArray(full(coil_subs)); Id_subs=gpuArray(Id_subs);
        
        % Run the calculation using backslash
        parfor n=1:numel(freq_grid)
            spectrum_subs=dot(coil_subs,((1i*L_subs+1i*freq_grid(n)*Id_subs)\rho0_subs));
            spectrum(n)=spectrum(n)+gather(spectrum_subs);
        end
        
    else
       
        % For large problems use GMRES
        if (size(rho0_subs,1)>5000)&&(size(rho0_subs,2)==1)
            
            % Inform the user
            report(spin_system,'using preconditioned CPU GMRES path...');
            
            % Get preconditioners
            [M1,M2]=ilu(1i*L_subs+1i*mean(freq_grid)*Id_subs,...
                        struct('type','crout','droptol',1e-3));
            report(spin_system,['nnz(L)='    num2str(nnz(L_subs)) ...
                                ', nnz(M1)=' num2str(nnz(M1))     ...
                                ', nnz(M2)=' num2str(nnz(M2))]);

            % MDCS diagnostics     
            parallel_profiler_start;                  
                          
            % Run the calculation in parallel
            parfor n=1:numel(freq_grid)
            
                % Run using preconditioned GMRES
                spectrum(n)=spectrum(n)+coil_subs'*gmres(1i*L_subs+1i*freq_grid(n)*Id_subs,...
                                        rho0_subs,10,1e-6,numel(rho0_subs),M1,M2);
                
            end
            
            % Get MDCS diagnostics       
            parallel_profiler_report;  
            
        else   
            
            % Inform the user
            report(spin_system,'using CPU backslash path...');
            
            % MDCS diagnostics     
            parallel_profiler_start; 
        
            % Run the calculation in parallel
            parfor n=1:numel(freq_grid)
                
                % Run using backslash on CPU
                spectrum(n)=spectrum(n)+dot(coil_subs,((1i*L_subs+1i*freq_grid(n)*Id_subs)\rho0_subs));
                
            end
            
            % Get MDCS diagnostics       
            parallel_profiler_report;
            
        end
    
    end
    
end

% Get the sampling rate implied by the frequency grid
sample_rate=abs(parameters.sweep(2)-parameters.sweep(1))*...
            parameters.npoints/(parameters.npoints-1);

% Match the unnormalised FFT amplitude convention
spectrum=spectrum*sample_rate;

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('spectral range must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||(numel(parameters.sweep)~=2)
    error('parameters.sweep vector must have two real elements');
end
if ~isfield(parameters,'npoints')
    error('number of points must be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<2)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be an integer greater than one.');
end
if ~isfield(parameters,'rho0')
    error('the initial state must be specified in parameters.rho0 variable.');
end
if (~isnumeric(parameters.rho0))||(~ismatrix(parameters.rho0))
    error('parameters.rho0 must be a numeric matrix.');
end
if size(parameters.rho0,1)~=size(H,2)
    error('parameters.rho0 must have the same number of rows as H.');
end
if ~isfield(parameters,'coil')
    error('the detection state must be specified in parameters.coil variable.');
end
if (~isnumeric(parameters.coil))||(~ismatrix(parameters.coil))
    error('parameters.coil must be a numeric matrix.');
end
if size(parameters.coil,1)~=size(H,2)
    error('parameters.coil must have the same number of rows as H.');
end
end

% "Our first rule here, Miss Taggart," he answered, "is that
%  one must always see for oneself."
%
% Ayn Rand, "Atlas Shrugged"

