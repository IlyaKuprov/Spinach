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
%       laxation matrix R should *not* be thermalized.
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
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    
    % Inform the user 
    report(spin_system,'projecting zeeman-hilb simulation into Liouville space...');
    
    % Project into Liouville space
    H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm');
    parameters.rho0=parameters.rho0(:); parameters.coil=parameters.coil(:);
    
    % Update the formalism setting
    spin_system.bas.formalism='zeeman-liouv';
    
end

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Compute subspace projectors
projectors=reduce(spin_system,L,parameters.coil);

% Loop over subspaces
for k=1:numel(projectors)
    
    % Project into the current subspace
    rho0_subs=projectors{k}'*parameters.rho0;
    coil_subs=projectors{k}'*parameters.coil;
    L_subs=projectors{k}'*L*projectors{k};
    Id_subs=speye(size(L_subs));
    
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
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('the initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('the detection state must be specified in parameters.coil variable.');
end
end

% "Our first rule here, Miss Taggart," he answered, "is that
%  one must always see for oneself."
%
% Ayn Rand, "Atlas Shrugged"

