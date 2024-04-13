% Computes resonance fields. For a Hamiltonian Hc+b*Hz, returns all 
% magnetic fields b for which the difference between two eigenvalues
% of Hc+b*Hz is equal to the frequency provided, and the transition
% moment across the specified operator Hmw is significant. Syntax:
%
%  [tf,tm,tw]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw)
%
% Parameters:
%
%     Iz  -  isotropic part of the laboratory frame Hamiltonian 
%            operator (Hilbert space) or commutation superopera-
%            tor (Liouville space, containing only Zeeman terms
%            at 1 Tesla
%
%     Qz  -  anisotropic part of the laboratory frame Hamiltoni-
%            an operator (Hilbert space) or commutation supero-
%            perator (Liouville space, containing only Zeeman 
%            terms at 1 Tesla
%
%     Ic  -  isotropic part of the laboratory frame Hamiltonian
%            operator (Hilbert space) or commutation superopera-
%            tor (Liouville space, containing all spin-spin cou-
%            plings, but no Zeeman terms
%
%     Qc  -  anisotropic part of the laboratory frame Hamiltoni-
%            an operator (Hilbert space) or commutation supero-
%            perator (Liouville space, containing all spin-spin
%            couplings, but no Zeeman terms
%
%     Hmw -  observable operator (Hilbert space) or observable
%            vector (Liouville space), without the amplitude 
%            prefactor
%
%     parameters.window   -  magnet field window, Tesla
%
%     parameters.mw_freq  -  microwave frequency, Hz
%
%     parameters.orientation - three Euler angles in radians
%                              specifying the system orientation
%
%     parameters.tm_tol   -  relative transition moment 
%                            tolerance
%
%     parameters.pp_tol   -  peak position tolerance in Tesla,
%                            this should be much smaller than
%                            the typical line width
%
%     parameters.rspt_order - perturbation theory order to use
%                             to account for the off-diagonal
%                             part of the Hamiltonian, Inf for
%                             exact diagonalisation
%
% Outputs:
%
%     tf     -  vector of transition fields in Tesla
%
%     tm     -  vector of transition moments
%
%     tw     -  vector of transition FWHMs in Tesla
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=eigenfields.m>

function [tf,tm,tw]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw)

% Check consistency
grumble(parameters,Iz,Qz,Ic,Qc);

% Recursive call for symmetry
if isfield(spin_system.bas,'irrep')
    
    % Preallocate irrep output blocks
    n_irreps=numel(spin_system.bas.irrep);
    tf=cell(n_irreps,1); tm=cell(n_irreps,1);
    tw=cell(n_irreps,1);
    
    % Loop over irreducible representations
    for n=1:n_irreps
        
        % Project Hamiltonian components
        P=spin_system.bas.irrep(n).projector;
        IzIrr=P'*Iz*P; IzIrr=(IzIrr+IzIrr')/2; 
        IcIrr=P'*Ic*P; IcIrr=(IcIrr+IcIrr')/2; 
        QzIrr=Qz; QcIrr=Qc;
        for k=1:numel(QzIrr)
            for m=1:numel(QzIrr{k})
                if ~isempty(QzIrr{k}{m})
                    QzIrr{k}{m}=P'*QzIrr{k}{m}*P;
                end
            end
        end
        for k=1:numel(QcIrr)
            for m=1:numel(QcIrr{k})
                if ~isempty(QcIrr{k}{m})
                    QcIrr{k}{m}=P'*QcIrr{k}{m}*P;
                end
            end
        end

        % Project the microwave operator
        switch spin_system.bas.formalism
            
            % Hilbert space
            case 'zeeman-hilb'
                
                % Hmw is a Hermitian matrix
                HmwIrr=P'*Hmw*P; HmwIrr=(HmwIrr+HmwIrr')/2;
                            
            % Liouville space
            case {'zeeman-liouv','sphten-liouv'}
                
                % Hmw is a vector
                HmwIrr=P'*Hmw;
                
            otherwise
                
                % Complain and bomb out
                error('unknown formalism specification.');
                
        end
           
        % Issue a recursive call
        spin_system_nosym=spin_system;
        spin_system_nosym.bas=rmfield(spin_system.bas,'irrep');
        [tf{n},tm{n},tw{n}]=eigenfields(spin_system_nosym,parameters,...
                                        IzIrr,QzIrr,IcIrr,QcIrr,HmwIrr);
                                    
    end
    
    % Concatenate irrep blocks
    tf=cell2mat(tf); tf=tf(:);
    tm=cell2mat(tm); tm=tm(:);
    tw=cell2mat(tw); tw=tw(:); return;
    
end

% Assemble and clean up Hamiltonians
Hc=Ic+orientation(Qc,parameters.orientation); 
Hz=Iz+orientation(Qz,parameters.orientation); 
Hc=(Hc+Hc')/2; Hz=(Hz+Hz')/2;

% Decide energy gap tolerance
engap_tol=abs(spin('E')*parameters.pp_tol);

% Microwave frequency
omega=2*pi*parameters.mw_freq;

% Pathway selection
switch spin_system.bas.formalism
    
    % Hilbert space formalism
    case 'zeeman-hilb'
        
        % Initial grid - trisect the window
        left_edge=min(parameters.window);
        right_edge=max(parameters.window);
        grid=[left_edge+0*(right_edge-left_edge)/3 ...
              left_edge+1*(right_edge-left_edge)/3 ...
              left_edge+2*(right_edge-left_edge)/3 ...
              left_edge+3*(right_edge-left_edge)/3];
          
        % Populate the initial grid
        E=zeros([size(Hz,1) numel(grid)]); dE=zeros(size(E)); 
        T=cell(1,numel(grid)); V=cell(1,numel(grid));
        for n=1:numel(grid)
            
            % Get sorted eigensystem
            [E(:,n),V{n}]=get_eigensys(parameters,grid(n)*Hz+Hc);
            
            % Hellmann-Feynman energy derivatives
            dE(:,n)=real(sum(conj(V{n}).*(Hz*V{n}),1));
            
            % Mark significant transition moments
            T{n}=abs(V{n}'*Hmw*V{n}).^2; T_max=max(T{n},[],'all');
            T{n}=sparse(T{n}>T_max*parameters.tm_tol);
            
        end
        
        % Paint all grid intervals dirty
        converged=false(numel(grid)-1,1);
        
        % Iterative grid trisection
        while any(~converged)
            
            % Start at the leftmost point of the grid
            new_E=E(:,1); new_dE=dE(:,1); new_T=T(1);
            new_grid=grid(1); new_conv=[]; new_V=V(1);
            
            % Inspect old grid intervals
            for n=2:numel(grid)
                
                % Check for prior convergence
                if converged(n-1)
                    
                    % Inherit the flag
                    new_conv(end+1)=true();  %#ok<AGROW>
                    
                else
                    
                    % Grid interval
                    dx=grid(n)-grid(n-1);
                    
                    % Predict and sort midpoint energies
                    EmP1=herm_spline(E(:,n-1),dE(:,n-1)*dx,...
                                     E(:,n),  dE(:,n)*dx, ...
                                     (1/3)*ones(size(E,1),1));
                    EmP2=herm_spline(E(:,n-1),dE(:,n-1)*dx,...
                                     E(:,n),  dE(:,n)*dx, ...
                                     (2/3)*ones(size(E,1),1));
                    EmP1=sort(EmP1); EmP2=sort(EmP2);
                                    
                    % Compute sorted midpoint energies
                    midp1=grid(n-1)+(1/3)*dx;
                    midp2=grid(n-1)+(2/3)*dx;
                    [Em1,V1]=get_eigensys(parameters,midp1*Hz+Hc);
                    [Em2,V2]=get_eigensys(parameters,midp2*Hz+Hc);
                    
                    % Hellmann-Feynman energy derivatives
                    dEm1=real(sum(conj(V1).*(Hz*V1),1));
                    dEm2=real(sum(conj(V2).*(Hz*V2),1));
                    
                    % Mark significant transition moments
                    Tm1=abs(V1'*Hmw*V1).^2; Tm1_max=max(Tm1,[],'all');
                    Tm2=abs(V2'*Hmw*V2).^2; Tm2_max=max(Tm2,[],'all');
                                        
                    % Add midpoints to the grid
                    new_grid(end+1)=midp1; new_grid(end+1)=midp2;  %#ok<AGROW>
                    new_dE(:,end+1)=dEm1;  new_dE(:,end+1)=dEm2;   %#ok<AGROW>
                    new_E(:,end+1)=Em1;    new_E(:,end+1)=Em2;     %#ok<AGROW>
                    new_T{end+1}=sparse(Tm1>Tm1_max*parameters.tm_tol);  %#ok<AGROW>
                    new_T{end+1}=sparse(Tm2>Tm2_max*parameters.tm_tol);  %#ok<AGROW>
                    new_V{end+1}=V1;  new_V{end+1}=V2;  %#ok<AGROW>
                    
                    % Check prediction accuracy
                    if (norm(EmP1-Em1,1)<engap_tol)&&...
                       (norm(EmP2-Em2,1)<engap_tol)
                        
                        % Paint new intervals clean
                        new_conv(end+1)=true();      %#ok<AGROW>
                        new_conv(end+1)=true();      %#ok<AGROW>
                        new_conv(end+1)=true();      %#ok<AGROW>
                        
                    else
                        
                        % Paint new intervals dirty
                        new_conv(end+1)=false();     %#ok<AGROW>
                        new_conv(end+1)=false();     %#ok<AGROW>
                        new_conv(end+1)=false();      %#ok<AGROW>
                         
                    end
                    
                end
                
                % Inherit the right point
                new_grid(end+1)=grid(n);   %#ok<AGROW>
                new_dE(:,end+1)=dE(:,n);   %#ok<AGROW>
                new_E(:,end+1)=E(:,n);     %#ok<AGROW>
                new_T{end+1}=T{n};         %#ok<AGROW>
                new_V{end+1}=V{n};         %#ok<AGROW>
                
            end
            
            % New grid becomes old
            E=new_E; dE=new_dE; grid=new_grid;
            T=new_T; V=new_V; converged=new_conv;
            
            % Complain and bomb out if the grid becomes too large
            if numel(grid)>1e3, error('grid construction failed'); end
            
        end
        
        % Get outputs started
        tf=[]; tm=[]; tw=[];
        
        % Loop over grid intervals
        for n=2:numel(grid)
            
            % Grid interval
            dx=grid(n)-grid(n-1);

            % Interval edge analysis
            interval_max=max(E(:,[n-1, n]),[],2);
            interval_min=min(E(:,[n-1, n]),[],2);

            % Screen by energy and transition moment
            deltaE_upper=interval_max'-interval_min;
            deltaE_lower=interval_min'-interval_max;
            promising_pairs=(omega<deltaE_upper)&...
                            (omega>deltaE_lower)&(T{n-1}|T{n});
            [source,destin]=find(promising_pairs);
            
            % Loop over source-destination pairs
            for k=1:numel(source)
                
                % Source and destination lines
                source_level=E(source(k),[n-1, n]);
                source_line=[source_level(2)-source_level(1); source_level(1)];
                destin_level=E(destin(k),[n-1, n]);
                destin_line=[destin_level(2)-destin_level(1); destin_level(1)];
                
                % Energies must be omega apart
                destin_line(2)=destin_line(2)-omega;

                % Solve the linear interpolant
                line_coeffs=destin_line-source_line;
                alpha=-line_coeffs(2)/line_coeffs(1);
                
                % Stay within the interval
                if (alpha>0)&&(alpha<1)

                    % Grid edge transition moments
                    tm_left=abs(V{n-1}(:,source(k))'*Hmw*V{n-1}(:,destin(k)))^2;
                    tm_right=abs(V{n}(:,source(k))'*Hmw*V{n}(:,destin(k)))^2;
                
                    % Store interpolated transition moment
                    tm(end+1)=(1-alpha)*tm_left+alpha*tm_right; %#ok<AGROW>
                  
                    % Store the eigenfield
                    tf(end+1)=grid(n-1)+alpha*dx; %#ok<AGROW>

                    % Get transition width (much to do here)
                    tw(end+1)=parameters.fwhm; %#ok<AGROW> 
                
                end
                
            end
        
        end
        
    % Liouville space formalism
    case {'zeeman-liouv','sphten-liouv'}
        
        % Get all transitions
        [uv,tf]=eig(omega*eye(size(Hc))-full(Hc),full(Hz),'vector');
        
        % Prune out instabilities
        hit_list=~isfinite(tf); uv(:,hit_list)=[]; tf(hit_list)=[];
        
        % Make sure tfs are real
        tf=real(tf);
        
        % Prune out-of-window transitions
        hit_list=(tf<min(parameters.window))|(tf>max(parameters.window));
        uv(:,hit_list)=[]; tf(hit_list)=[];

        % Normalise dyadics
        uv=uv./sqrt(diag(uv'*uv)');
        
        % Compute transition moments
        tm=abs(Hmw'*uv).^2;

        % Get transition widths (much to do here)
        tw=ones(size(tm))*parameters.fwhm;

    otherwise
        
        % Complain and bomb out
        error('unexpected formalism specification');
    
end

% Prune insignificant transition moments
tm_max=max(tm,[],'all');
hit_list=(tm<tm_max*parameters.tm_tol);
tf(hit_list)=[]; tm(hit_list)=[]; tw(hit_list)=[];

% Reshape into columns and sort
[tf,idx]=sort(tf(:)); 
tm=tm(:); tm=tm(idx);
tw=tw(:); tw=tw(idx);

end

% Eigensystem with an RSPT switch
function [E,V]=get_eigensys(parameters,H)

% Decide the method
switch parameters.rspt_order

    case {1,2,3,4}

        % Diag and off-diag Hamiltonian
        H_diag=diag(H); H_offd=(H+H')/2; 
        H_offd(logical(speye(size(H))))=0;

        % Use perturbation theory
        [E,V]=rspert(H_diag,H_offd,parameters.rspt_order);

    case Inf

        % Full diagonalisation
        [V,E]=eig(full(H),'vector');

    otherwise

        % Complain and bomb out
        error('unsupported perturbation theory order.');

end

% Eigensystem sorting
[E,idx]=sort(E); V=V(:,idx);

end

% Consistency enforcement
function grumble(parameters,Iz,Qz,Ic,Qc)
if (~isnumeric(Ic))||(size(Ic,1)~=size(Ic,2))||...
   (~isnumeric(Iz))||(size(Iz,1)~=size(Iz,2))||...
   (size(Ic,1)~=size(Iz,1))
    error('Ic and Iz must be square matrices of the same size.');
end
if (~ishermitian(Ic))||(~ishermitian(Iz))
     error('Ic and Iz must be Hermitian.');
end
if (~iscell(Qz))||(~iscell(Qc))
    error('Qz and Qc must be cell arrays produced by hamiltonian() function.');
end
if ~isfield(parameters,'mw_freq')
    error('resonance frequency must be supplied in parameters.mw_freq field.');
end
if ~isfield(parameters,'window')
    error('field window must be supplied in parameters.window field.');
end
if ~isfield(parameters,'pp_tol')
    error('peak position tolerance must be supplied in parameters.pp_tol field.');
end
if ~isfield(parameters,'tm_tol')
    error('transition moment tolerance must be supplied in parameters.tm_tol field.');
end
if ~isfield(parameters,'orientation')
    error('system orientation must be supplied in parameters.orientation field.');
end
if (~isnumeric(parameters.mw_freq))||...
   (~isreal(parameters.mw_freq))||...
   (~isscalar(parameters.mw_freq))
    error('parameters.mw_freq must be a real scalar.');
end
if (~isnumeric(parameters.pp_tol))||...
   (~isreal(parameters.pp_tol))||...
   (~isscalar(parameters.pp_tol))
    error('parameters.pp_tol must be a real scalar.');
end
if (~isnumeric(parameters.tm_tol))||...
   (~isreal(parameters.tm_tol))||...
   (~isscalar(parameters.tm_tol))
    error('parameters.tm_tol must be a real scalar.');
end
if (~isnumeric(parameters.window))||...
   (~isreal(parameters.window))||...
   (numel(parameters.window)~=2)
    error('parameters.window must have two real elements.');
end
if (~isnumeric(parameters.orientation))||...
   (~isreal(parameters.orientation))||...
   (numel(parameters.orientation)~=3)
    error('parameters.orientation must have three real elements.');
end
end

% But even in my life I saw the leaching of spirit. A surfeit of 
% honey cloys the tongue; a surfeit of wine addles the brain; so
% a surfeit of ease guts a man of strength. Light, warmth, food,
% water were free to all men, and gained by a minimum of effort.
% So the people of Ampridatvir, released from toil, gave increa-
% sing attention to faddishness, perversity, and the occult.
%
% Jack Vance, "The Dying Earth"

