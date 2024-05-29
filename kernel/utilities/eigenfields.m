% Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all 
% magnetic fields B for which the difference between two eigenvalues
% of Hc+B*Hz is equal to the frequency provided, and the transition
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
%     pd     -  vector of energy level population differences
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=eigenfields.m>

function [tf,tm,tw,pd]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw)

% Check consistency
grumble(parameters,Iz,Qz,Ic,Qc);

% Get frequency gap tolerance in rad/s
frq_gap_tol=abs(spin('E')*parameters.pp_tol);

% Assemble Zeeman and coupling Hamiltonians
Hz=Iz+orientation(Qz,parameters.orientation); 
Hc=Ic+orientation(Qc,parameters.orientation); 

% Get MW frequency in rad/s
omega=2*pi*parameters.mw_freq;

% Pathway selection
switch spin_system.bas.formalism
    
    % Hilbert space
    case 'zeeman-hilb'
        
        % Initial field grid: trisect the window
        left_edge=min(parameters.window);
        right_edge=max(parameters.window);
        grid=[left_edge+0*(right_edge-left_edge)/3 ...
              left_edge+1*(right_edge-left_edge)/3 ...
              left_edge+2*(right_edge-left_edge)/3 ...
              left_edge+3*(right_edge-left_edge)/3];
          
        % Set up field grid data structures
        E= zeros([size(Hz,1) numel(grid)]); % Level energies
        LP=zeros([size(Hz,1) numel(grid)]); % Level populations
        dE=zeros([size(Hz,1) numel(grid)]); % dE/dB derivatives
        T=cell(1,numel(grid));              % Transition moments

        % Populate initial grid
        for n=1:numel(grid)

            % Sorted (ascending) energies, dE/dB derivatives, transition moments, level pops
            [E(:,n),~,dE(:,n),T{n},LP(:,n)]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,grid(n));
            
        end
        
        % Paint all grid intervals dirty
        converged=false(numel(grid)-1,1);
        
        % Iterative grid trisection
        while any(~converged)
            
            % Start at the leftmost point of the grid
            new_E=E(:,1); new_dE=dE(:,1); new_T=T(1);
            new_grid=grid(1); new_conv=[]; new_LP=LP(:,1);
            
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
                                    
                    % Get eigensystem information for the grid midpoints
                    midp1=grid(n-1)+(1/3)*dx; midp2=grid(n-1)+(2/3)*dx;
                    [Em1,~,dEm1,Tm1,LPm1]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,midp1);
                    [Em2,~,dEm2,Tm2,LPm2]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,midp2);

                    % Append midpoints to the new grid
                    new_grid((end+1):(end+2))=[midp1 midp2];
                    new_dE(:,(end+1):(end+2))=[dEm1 dEm2];
                    new_LP(:,(end+1):(end+2))=[LPm1 LPm2];
                    new_E(:,(end+1):(end+2))=[Em1 Em2];
                    new_T((end+1):(end+2))={Tm1 Tm2};

                    % Check prediction accuracy
                    if (norm(EmP1-Em1,1)<frq_gap_tol)&&...
                       (norm(EmP2-Em2,1)<frq_gap_tol)
                        
                        % Paint new intervals clean
                        new_conv((end+1):(end+3))=true();
                        
                    else
                        
                        % Paint new intervals dirty
                        new_conv((end+1):(end+3))=false();     
                         
                    end
                    
                end
                
                % Inherit right point from old grid
                new_grid(end+1)=grid(n);   %#ok<AGROW>
                new_dE(:,end+1)=dE(:,n);   %#ok<AGROW>
                new_E(:,end+1)=E(:,n);     %#ok<AGROW>
                new_T{end+1}=T{n};         %#ok<AGROW>
                new_LP(:,end+1)=LP(:,n);   %#ok<AGROW>
                
            end
            
            % New grid becomes old
            E=new_E; dE=new_dE; grid=new_grid;
            T=new_T; LP=new_LP; converged=new_conv;

            % Complain and bomb out if the field grid becomes too large
            if numel(grid)>1e3, error('field grid construction failed'); end
            
        end
        
        % Get outputs started
        tf=[]; tm=[]; tw=[]; pd=[];

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
            promising_pairs=(omega<deltaE_upper)&(omega>deltaE_lower)&...
                            (T{n-1}>parameters.tm_tol|T{n}>parameters.tm_tol);
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

                % Use linear interpolation first
                line_coeffs=destin_line-source_line;
                alpha=-line_coeffs(2)/line_coeffs(1);

                % Check level crossing location
                if (alpha<-0.1)||(alpha>1.1)
                    
                    % Outside
                    continue;

                else

                    % Get spline coefficients (x^3 -> x^0) for energies in this interval
                    source_spline=[dE(source(k),n-1)*dx E(source(k),n-1) ...
                                   dE(source(k),n)*dx   E(source(k),n)]*[ 1 -2  1  0; 
                                                                          2 -3  0  1; 
                                                                          1 -1  0  0; 
                                                                         -2  3  0  0];
                    destin_spline=[dE(destin(k),n-1)*dx E(destin(k),n-1) ...
                                   dE(destin(k),n)*dx   E(destin(k),n)]*[ 1 -2  1  0; 
                                                                          2 -3  0  1; 
                                                                          1 -1  0  0; 
                                                                         -2  3  0  0];
                    % Energies must be omega apart
                    destin_spline(4)=destin_spline(4)-omega;

                    % Get and stabilise cubic equation coefficients
                    cubic_coeffs=destin_spline-source_spline;
                    cubic_coeffs=cubic_coeffs/max(abs(cubic_coeffs));

                    % Differentiate the spline and store coefficients as (x^2 -> x^0)
                    deriv_coeffs=[3*cubic_coeffs(1) 2*cubic_coeffs(2) 1*cubic_coeffs(3)];

                    % One iteration of Newton's method to refine the root
                    numer=cubic_coeffs(1)*alpha^3+cubic_coeffs(2)*alpha^2+...
                          cubic_coeffs(3)*alpha  +cubic_coeffs(4);
                    denom=deriv_coeffs(1)*alpha^2+deriv_coeffs(2)*alpha  +...
                          deriv_coeffs(3);
                    alpha=alpha-numer/denom;

                end
                    
                % Check level crossing location
                if (alpha<0)||(alpha>1)
                    
                    % Outside
                    continue;

                else

                    % Interval edge transition moments
                    tm_left=T{n-1}(source(k),destin(k));
                    tm_right=T{n}(source(k),destin(k));

                    % Interval edge population differences
                    pd_left=LP(source(k),n-1)-LP(destin(k),n-1);
                    pd_right=LP(source(k),n)-LP(destin(k),n);
                
                    % Store interpolated transition moment
                    tm(end+1)=(1-alpha)*tm_left+alpha*tm_right; %#ok<AGROW>

                    % Store interpolated population difference
                    pd(end+1)=(1-alpha)*pd_left+alpha*pd_right; %#ok<AGROW>
                  
                    % Store the transition field
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
        
        % Prune out unphysical results
        hit_list=~isfinite(tf); tf(hit_list)=[];  
        uv(:,hit_list)=[]; tf=real(tf);
        
        % Prune out-of-window transitions
        hit_list=(tf<min(parameters.window))|...
                 (tf>max(parameters.window));
        uv(:,hit_list)=[]; tf(hit_list)=[];

        % Normalise dyadics
        uv=uv./sqrt(diag(uv'*uv)');
        
        % Compute transition moments
        tm=abs(Hmw'*uv).^2;

        % Unit pop diffs for now
        pd=ones(size(tm));

        % Get transition widths (much to do here)
        tw=ones(size(tm))*parameters.fwhm;

    otherwise
        
        % Complain and bomb out
        error('unexpected formalism specification');
    
end

% Prune insignificant transition moments
hit_list=(tm<parameters.tm_tol);
tf(hit_list)=[]; tm(hit_list)=[]; 
tw(hit_list)=[]; pd(hit_list)=[];

% Reshape into columns and sort
[tf,idx]=sort(tf(:)); 
tm=tm(:); tm=tm(idx);
tw=tw(:); tw=tw(idx);
pd=pd(:); pd=pd(idx);

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

