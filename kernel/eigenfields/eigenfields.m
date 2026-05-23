% Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all 
% magnetic fields B for which the difference between two eigenvalues
% of Hc+B*Hz is equal to the frequency provided, and the transition
% moment across the specified operator Hmw is significant. Syntax:
%
%  [tf,tm,tw,pd,ti,tj]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw)
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
%     parameters.fwhm     -  transition full width at half
%                            maximum, Tesla
%
%     parameters.rspt_order - perturbation theory order to use
%                             to account for the off-diagonal
%                             part of the Hamiltonian, Inf for
%                             exact diagonalisation
%
% Outputs:
%
%     tran.tf     -  vector of transition fields in Tesla
%
%     tran.tm     -  vector of transition moments
%
%     tran.tw     -  vector of transition FWHMs in Tesla
%
%     tran.pd     -  vector of energy level population differences
%
%     tran.ti     -  transition identity array, one row per transition
%
%     tran.tj     -  vector of scaled field-sweep Jacobians
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=eigenfields.m>

function tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)

% Check consistency
grumble(spin_system,parameters,Hz,Hc);

% Get frequency gap tolerance in rad/s
frq_gap_tol=abs(spin('E')*parameters.pp_tol);

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
        V=cell(1,numel(grid));              % Level eigenvectors
        T=cell(1,numel(grid));              % Transition moments

        % Populate field grid
        for n=1:numel(grid)

            % Sorted (ascending) energies, dE/dB derivatives, transition moments, level pops
            [E(:,n),V{n},dE(:,n),T{n},LP(:,n)]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,grid(n));
            
        end
        
        % Paint all intervals dirty
        converged=false(numel(grid)-1,1);
        
        % Iterative trisection
        while any(~converged)
            
            % Start at the leftmost point of the field grid
            new_grid=grid(1); new_E=E(:,1); new_dE=dE(:,1);
            new_T=T(1); new_V=V(1); new_LP=LP(:,1); new_conv=[];
            
            % Inspect old intervals
            for n=2:numel(grid)
                
                % Check for prior convergence
                if converged(n-1)
                    
                    % Inherit the flag
                    new_conv(end+1)=true();  %#ok<AGROW>
                    
                else
                    
                    % Field grid interval
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
                    [Em1,Vm1,dEm1,Tm1,LPm1]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,midp1);
                    [Em2,Vm2,dEm2,Tm2,LPm2]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,midp2);

                    % Append midpoints to the new field grid
                    new_grid((end+1):(end+2))=[midp1 midp2];
                    new_dE(:,(end+1):(end+2))=[dEm1 dEm2];
                    new_LP(:,(end+1):(end+2))=[LPm1 LPm2];
                    new_E(:,(end+1):(end+2))=[Em1 Em2];
                    new_T((end+1):(end+2))={Tm1 Tm2};
                    new_V((end+1):(end+2))={Vm1 Vm2};

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
                new_V{end+1}=V{n};         %#ok<AGROW>
                new_LP(:,end+1)=LP(:,n);   %#ok<AGROW>
                
            end
            
            % New field grid becomes old
            E=new_E; dE=new_dE; grid=new_grid;
            T=new_T; V=new_V; LP=new_LP; converged=new_conv;

            % Complain and bomb out if the field grid becomes too large
            if numel(grid)>1e3, error('field grid construction failed'); end
            
        end
        
        % Reorder eigenstates by maximum overlap between field knots
        for n=2:numel(grid)

            % Build a one-to-one state assignment
            state_ovlp=abs(V{n-1}'*V{n});
            perm=zeros(size(E,1),1);
            state_hit=false(size(E,1),1);
            for k=1:size(E,1)

                % Match the largest remaining eigenvector overlap
                [max_ovlp,linear_idx]=max(state_ovlp(:));
                if max_ovlp<sqrt(eps), error('state tracking failed'); end
                [left_idx,right_idx]=ind2sub(size(state_ovlp),linear_idx);
                perm(left_idx)=right_idx;
                state_hit(left_idx)=true;
                state_ovlp(left_idx,:)=-Inf;
                state_ovlp(:,right_idx)=-Inf;

            end
            if any(~state_hit), error('state tracking failed'); end

            % Reorder right-knot quantities into left-knot labelling
            E(:,n)=E(perm,n); dE(:,n)=dE(perm,n);
            LP(:,n)=LP(perm,n); V{n}=V{n}(:,perm);
            T{n}=T{n}(perm,perm);

        end

        % Get the relative transition moment cut-off
        tm_scale=0;
        for n=1:numel(T)
            tm_scale=max([tm_scale; T{n}(:)]);
        end
        tm_cutoff=parameters.tm_tol*tm_scale;

        % Mark active level pairs
        pair_active=false(size(E,1));
        for n=1:numel(T)
            pair_active=pair_active|(T{n}>tm_cutoff);
        end

        % Compile the active level pair list
        [source,destin]=find(pair_active);

        % Get the outputs started
        tran.tf=[]; tran.tm=[]; tran.tw=[]; 
        tran.pd=[]; tran.tj=[]; tran.ti=zeros(0,3);

        % Initialise branch counters
        pair_branch=zeros(size(E,1));

        % Loop over grid intervals
        for n=2:numel(grid)
            
            % Grid interval
            dx=grid(n)-grid(n-1);

            % State overlap between field knots
            if isempty(source), continue; end
            state_ovlp=abs(V{n-1}'*V{n}).^2;
            
            % Loop over source-destination pairs
            for k=1:numel(source)
                
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

                % Get the cubic equation coefficients
                gap_spline=destin_spline-source_spline;
                poly_coeffs=gap_spline;

                % Scale for numerical stability
                poly_scale=max(abs(poly_coeffs));
                if poly_scale==0, continue; end
                poly_coeffs=poly_coeffs/poly_scale;

                % Drop leading numerical zeros
                lead_idx=find(abs(poly_coeffs)>sqrt(eps),1,'first');
                if isempty(lead_idx), continue; end
                poly_coeffs=poly_coeffs(lead_idx:end);

                % Find all real roots inside the field interval
                root_list=roots(poly_coeffs); root_tol=sqrt(eps);
                root_list=root_list(abs(imag(root_list))<root_tol);
                root_list=real(root_list);
                root_list=root_list((root_list>=-root_tol)&...
                                    (root_list<=1+root_tol));
                root_list=min(max(root_list,0),1);

                % Add tangent roots at spline extrema
                turn_list=roots(polyder(gap_spline));
                turn_list=turn_list(abs(imag(turn_list))<root_tol);
                turn_list=real(turn_list);
                turn_list=turn_list((turn_list>=-root_tol)&...
                                    (turn_list<=1+root_tol));
                turn_list=min(max(turn_list,0),1);
                turn_list=turn_list(abs(polyval(gap_spline,turn_list))<frq_gap_tol);
                root_list=sort([root_list(:).' turn_list(:).']);
                if isempty(root_list), continue; end
                root_list=root_list([true abs(diff(root_list))>root_tol]);

                % Check state labelling stability
                stable_states=(state_ovlp(source(k),source(k))>0.5)&&...
                              (state_ovlp(destin(k),destin(k))>0.5);

                % Loop over all roots in this field interval
                for q=1:numel(root_list)

                    % Get the root coordinate
                    alpha=root_list(q);
                    if (alpha<=root_tol)&&(n>2), continue; end

                    % Interval edge transition moments
                    tm_left=T{n-1}(source(k),destin(k));
                    tm_right=T{n}(source(k),destin(k));

                    % Interval edge population differences
                    pd_left=LP(source(k),n-1)-LP(destin(k),n-1);
                    pd_right=LP(source(k),n)-LP(destin(k),n);

                    % Interpolate transition properties
                    tm_base=(1-alpha)*tm_left+alpha*tm_right;
                    pd_base=(1-alpha)*pd_left+alpha*pd_right;

                    % Get the transition frequency slope
                    jac_slope=polyval(polyder(gap_spline),alpha)/dx;

                    % Rediagonalise at unstable roots
                    if ~stable_states
                        reson_field=grid(n-1)+alpha*dx;
                        [~,V_root,dE_root,T_root,LP_root]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,reson_field);

                        % Match root eigenvectors to the left knot
                        root_ovlp=abs(V{n-1}'*V_root);
                        perm=zeros(size(E,1),1);
                        state_hit=false(size(E,1),1);
                        for m=1:size(E,1)

                            % Match the largest remaining eigenvector overlap
                            [max_ovlp,linear_idx]=max(root_ovlp(:));
                            if max_ovlp<sqrt(eps), error('state tracking failed'); end
                            [left_idx,right_idx]=ind2sub(size(root_ovlp),linear_idx);
                            perm(left_idx)=right_idx;
                            state_hit(left_idx)=true;
                            root_ovlp(left_idx,:)=-Inf;
                            root_ovlp(:,right_idx)=-Inf;

                        end
                        if any(~state_hit), error('state tracking failed'); end

                        % Reorder root quantities into left-knot labelling
                        dE_root=dE_root(perm);
                        T_root=T_root(perm,perm);
                        LP_root=LP_root(perm);
                        tm_base=T_root(source(k),destin(k));
                        pd_base=LP_root(source(k))-LP_root(destin(k));
                        jac_slope=dE_root(destin(k))-dE_root(source(k));
                    end

                    % Scale the field-sweep Jacobian to electron EPR order
                    if jac_slope==0
                        tj_base=Inf;
                    else
                        tj_base=abs(spin('E'))/abs(jac_slope);
                    end

                    % Update the branch count for this level pair
                    pair_branch(source(k),destin(k))=pair_branch(source(k),destin(k))+1;

                    % Store interpolated transition moment
                    tran.tm(end+1)=tm_base; 

                    % Store population difference
                    tran.pd(end+1)=pd_base; 

                    % Store scaled field-sweep Jacobian
                    tran.tj(end+1)=tj_base; 

                    % Store the transition field
                    tran.tf(end+1)=grid(n-1)+alpha*dx;

                    % Get transition width (much to do here)
                    tran.tw(end+1)=parameters.fwhm; 

                    % Store transition identity
                    tran.ti(end+1,:)=[source(k) destin(k) pair_branch(source(k),destin(k))];

                end
                
            end
        
        end
        
    % Liouville space formalism
    case {'zeeman-liouv','sphten-liouv'}
        
        % Set up the generalised eigenfield pencil
        left_mat=omega*eye(size(Hc))-full(Hc);
        right_mat=full(Hz);

        % Get all transitions
        [uv,tf]=eig(left_mat,right_mat,'vector');
        tran.tf=tf(:).';
        
        % Prune out unphysical results
        uv_norm=sqrt(sum(abs(uv).^2,1));
        resid=vecnorm(left_mat*uv-right_mat*(uv.*tf.'),2,1);
        resid_scale=norm(left_mat,2)+norm(right_mat,2)*max(1,abs(tf(:).'));
        complex_tf=abs(imag(tf(:).'))>1e-8*max(1,abs(real(tf(:).')));
        hit_list=(~isfinite(real(tf(:).')))|...
                 (~isfinite(imag(tf(:).')))|complex_tf|...
                 (~isfinite(uv_norm))|(uv_norm<sqrt(eps))|...
                 (resid>1e-8*resid_scale);
        tran.tf(hit_list)=[]; uv(:,hit_list)=[];
        tran.tf=real(tran.tf);
        
        % Prune out-of-window transitions
        hit_list=(tran.tf<min(parameters.window))|...
                 (tran.tf>max(parameters.window));
        uv(:,hit_list)=[]; tran.tf(hit_list)=[];

        % Normalise dyadics
        uv_norm=sqrt(sum(abs(uv).^2,1)); uv=uv./uv_norm;
        
        % Compute transition moments
        tran.tm=abs(Hmw'*uv).^2;

        % Compute scaled field-sweep Jacobians
        jac_slope=real(sum(conj(uv).*(right_mat*uv),1));
        tran.tj=abs(spin('E'))./abs(jac_slope);

        % Unit pop diffs for now
        tran.pd=ones(size(tran.tm));

        % Get transition widths (much to do here)
        tran.tw=ones(size(tran.tm))*parameters.fwhm;

        % Use the generalised eigenvector ordinal as transition identity
        tran.ti=(1:numel(tran.tf)).';

    otherwise
        
        % Complain and bomb out
        error('unexpected formalism specification');
    
end

% Prune insignificant transition moments
if isempty(tran.tm)
    hit_list=false(size(tran.tm));
elseif exist('tm_cutoff','var')
    hit_list=(tran.tm<tm_cutoff);
else
    tm_cutoff=parameters.tm_tol*max(tran.tm(:));
    hit_list=(tran.tm<tm_cutoff);
end
tran.tf(hit_list)=[]; tran.tm(hit_list)=[];
tran.tw(hit_list)=[]; tran.pd(hit_list)=[]; 
tran.tj(hit_list)=[]; tran.ti(hit_list,:)=[];

% Reshape into columns and sort
[tran.tf,idx]=sort(tran.tf(:)); 
tran.tm=tran.tm(:); tran.tm=tran.tm(idx);
tran.tw=tran.tw(:); tran.tw=tran.tw(idx);
tran.pd=tran.pd(:); tran.pd=tran.pd(idx);
tran.tj=tran.tj(:); tran.tj=tran.tj(idx);
tran.ti=tran.ti(idx,:);

% Blank coordinates for later
tran.xyz=[NaN; NaN; NaN];

end

% Consistency enforcement
function grumble(spin_system,parameters,Hz,Hc)
if (~isnumeric(Hc))||(size(Hc,1)~=size(Hc,2))||...
   (~isnumeric(Hz))||(size(Hz,1)~=size(Hz,2))||...
   (size(Hc,1)~=size(Hz,1))
    error('Hc and Hz must be square matrices of the same size.');
end
if (~ishermitian(Hc))||(~ishermitian(Hz))
     error('Hc and Hz must be Hermitian.');
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
if ~isfield(parameters,'fwhm')
    error('transition FWHM must be supplied in parameters.fwhm field.');
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
if (~isnumeric(parameters.fwhm))||...
   (~isreal(parameters.fwhm))||...
   (~isscalar(parameters.fwhm))||...
   (parameters.fwhm<=0)
    error('parameters.fwhm must be a positive real scalar.');
end
if (~isnumeric(parameters.window))||...
   (~isreal(parameters.window))||...
   (numel(parameters.window)~=2)
    error('parameters.window must have two real elements.');
end
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    if ~isfield(parameters,'rspt_order')
        error('perturbation theory order must be supplied in parameters.rspt_order field.');
    end
    if (~isnumeric(parameters.rspt_order))||(~isreal(parameters.rspt_order))||...
       (~isscalar(parameters.rspt_order))||((mod(parameters.rspt_order,1)~=0)&&...
       (~isinf(parameters.rspt_order)))||(parameters.rspt_order<0)
        error('parameters.rspt_order must be a non-negative integer or Inf.');
    end
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
