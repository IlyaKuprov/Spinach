% Relaxation superoperator. Syntax:
%
%                R=relaxation(spin_system,euler_angles)
%
% Parameters:
%
%     euler_angles  - if this parameter is skipped, the isotropic
%                     part of the relaxation superoperator is re-
%                     turned; when this parameter is specified,
%                     those theories that support relaxation ani-
%                     sotropy start taking spin system orientati-
%                     on into account
%                    
% Outputs:
%
%     R             - relaxation superoperator
%
% Note: a variety of relaxation theories are supported, see the relax-
%       ation theory parameters section of the online manual.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk 
%
% <https://spindynamics.org/wiki/index.php?title=relaxation.m>

function R=relaxation(spin_system,euler_angles)

% Set defaults
spin_system=defaults(spin_system);

% Check consistency
grumble(spin_system);

% Get the matrix going
R=mprealloc(spin_system,0);

% Do nothing if none specified
if isempty(spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'relaxation superoperator set to zero.');
    
    % Exit the function
    return;
    
end

% Add extended T1/T2 model terms
if ismember('t1_t2',spin_system.rlx.theories)
    
    % Call the extended T1/T2 model function
    if exist('euler_angles','var')
        report(spin_system,'adding anisotropic extended T1/T2 model terms...');
        [R1,R2]=rlx_t1_t2(spin_system,euler_angles);
    else
        report(spin_system,'adding isotropic extended T1/T2 model terms...');
        [R1,R2]=rlx_t1_t2(spin_system); 
    end
    R=R+R1+R2;
    
end

% Add Redfield terms
if ismember('redfield',spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'Redfield theory with user-supplied correlation times.');

    % Catch zero correlation times
    for n=1:numel(spin_system.rlx.tau_c)
        if ~any(spin_system.rlx.tau_c{n},'all')
            error('zero rotational correlation times are not supported.');
        end
    end

    % Get the rotational basis, including the non-secular terms
    report(spin_system,'computing the lab frame Hamiltonian superoperator...');
    [L0,Q]=hamiltonian(assume(spin_system,'labframe')); %#ok<ASGLU>
    
    % Compute Redfield integral
    if isworkernode||ismember('asyredf',spin_system.sys.disable)
        report(spin_system,'serial evaluation path...');
        redfield_integral_serial;
    else
        report(spin_system,'parallel evaluation path...');
        redfield_integral_async;
    end
    
    % Catch abuse of Redfield theory
    max_rate=max(abs(diag(R)));
    for n=1:numel(spin_system.rlx.tau_c)
        if (1/max_rate)<spin_system.rlx.tau_c{n}
            report(spin_system,['1/max(diag(R)) = ' num2str(1/max_rate) ...
                                ', tau_c = ' num2str(spin_system.rlx.tau_c{n})]);
            error('T1,2>>tau_c validity condition violation in Redfield theory');
        end
    end
    
end
        
% Add Lindblad terms
if ismember('lindblad',spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'Lindblad theory with user-supplied R1 and R2 rates.');
    
    % Loop over spins
    for n=1:spin_system.comp.nspins
        
        % Get the basic superoperators
        Lp_left=operator(spin_system,{'L+'},{n},'left');
        Lp_right=operator(spin_system,{'L+'},{n},'right');
        Lm_left=operator(spin_system,{'L-'},{n},'left');
        Lm_right=operator(spin_system,{'L-'},{n},'right');
        Lz_left=operator(spin_system,{'Lz'},{n},'left');
        Lz_right=operator(spin_system,{'Lz'},{n},'right');
        
        % Get the rates
        Rpm=spin_system.rlx.lind_r1_rates(n)/4;
        Rmp=spin_system.rlx.lind_r1_rates(n)/4;
        Rzz=spin_system.rlx.lind_r2_rates(n)-...
            spin_system.rlx.lind_r1_rates(n)/2;
        
        % Build the superoperator
        R=R+Rpm*(2*Lp_left*Lm_right-Lm_left*Lp_left-Lp_right*Lm_right)+...
            Rmp*(2*Lm_left*Lp_right-Lp_left*Lm_left-Lm_right*Lp_right)+...
            Rzz*(2*Lz_left*Lz_right-Lz_left*Lz_left-Lz_right*Lz_right);
        
    end
    
    % Inform the user
    report(spin_system,'Lindblad theory terms have been incorporated.');
    
end
    
% Add Weizmann DNP theory terms
if ismember('weizmann',spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'Weizmann DNP theory with user-supplied R1e, R2e, R1n, R2n, R1d and R2d rates.');
    
    % Get electron superoperators
    Ep_L=operator(spin_system,'L+','electrons','left');
    Ep_R=operator(spin_system,'L+','electrons','right');
    Em_L=operator(spin_system,'L-','electrons','left');
    Em_R=operator(spin_system,'L-','electrons','right');
    Ez_L=operator(spin_system,'Lz','electrons','left');
    Ez_R=operator(spin_system,'Lz','electrons','right');
    
    % Get electron rates
    Rpm=spin_system.rlx.weiz_r1e/4; 
    Rmp=spin_system.rlx.weiz_r1e/4;
    Rzz=spin_system.rlx.weiz_r2e-...
        spin_system.rlx.weiz_r1e/2;
    
    % Add Lindblad type terms
    R=R+Rpm*(2*Ep_L*Em_R-Em_L*Ep_L-Ep_R*Em_R)+...
        Rmp*(2*Em_L*Ep_R-Ep_L*Em_L-Em_R*Ep_R)+...
        Rzz*(2*Ez_L*Ez_R-Ez_L*Ez_L-Ez_R*Ez_R);
    
    % Get nuclear superoperators
    Np_L=operator(spin_system,'L+','nuclei','left');
    Np_R=operator(spin_system,'L+','nuclei','right');
    Nm_L=operator(spin_system,'L-','nuclei','left');
    Nm_R=operator(spin_system,'L-','nuclei','right');
    Nz_L=operator(spin_system,'Lz','nuclei','left');
    Nz_R=operator(spin_system,'Lz','nuclei','right');
    
    % Get nuclear rates
    Rpm=spin_system.rlx.weiz_r1n/4;
    Rmp=spin_system.rlx.weiz_r1n/4;
    Rzz=spin_system.rlx.weiz_r2n-...
        spin_system.rlx.weiz_r1n/2;
    
    % Add Lindblad type terms
    R=R+Rpm*(2*Np_L*Nm_R-Nm_L*Np_L-Np_R*Nm_R)+...
        Rmp*(2*Nm_L*Np_R-Np_L*Nm_L-Nm_R*Np_R)+...
        Rzz*(2*Nz_L*Nz_R-Nz_L*Nz_L-Nz_R*Nz_R);
    
    % Add dipolar cross-relaxation terms
    for n=1:spin_system.comp.nspins
        for k=1:spin_system.comp.nspins
            if n~=k
            
                % Process longitudinal terms
                if spin_system.rlx.weiz_r1d(n,k)~=0
                
                    % Generate flip-flop superoperator
                    FF_L=operator(spin_system,{'L+','L-'},{n,k},'left')+...
                         operator(spin_system,{'L-','L+'},{n,k},'left');
                    FF_R=operator(spin_system,{'L+','L-'},{n,k},'right')+...
                         operator(spin_system,{'L-','L+'},{n,k},'right');
                
                    % Add Lindblad type term
                    R=R+spin_system.rlx.weiz_r1d(n,k)*(2*FF_L*FF_R-FF_L*FF_L-FF_R*FF_R)/2;
                
                end
            
                % Process transverse terms
                if spin_system.rlx.weiz_r2d(n,k)~=0
                
                    % Generate ZZ superoperator
                    ZZ_L=operator(spin_system,{'Lz','Lz'},{n,k},'left');
                    ZZ_R=operator(spin_system,{'Lz','Lz'},{n,k},'right');
                
                    % Add Lindblad type term
                    R=R+spin_system.rlx.weiz_r2d(n,k)*(2*ZZ_L*ZZ_R-ZZ_L*ZZ_L-ZZ_R*ZZ_R)/2;
                
                end

            end
        end
    end
    
    % Inform the user
    report(spin_system,'Weizmann DNP relaxation theory terms have been incorporated.');
    
end
        
% Add Nottingham DNP theory terms
if ismember('nottingham',spin_system.rlx.theories)
        
    % Update the user
    report(spin_system,'Nottingham DNP theory with user-supplied R1e, R1n, R2e and R2n rates.');
    
    % Set shorthands
    r1e=spin_system.rlx.nott_r1e; r2e=spin_system.rlx.nott_r2e;
    r1n=spin_system.rlx.nott_r1n; r2n=spin_system.rlx.nott_r2n;
    
    % Find electrons
    electrons=find(strcmp('E',spin_system.comp.isotopes));
    
    % Get basic electron operators
    S1p=operator(spin_system,{'L+'},{electrons(1)});
    S1z=operator(spin_system,{'Lz'},{electrons(1)});
    S2p=operator(spin_system,{'L+'},{electrons(2)});
    S2z=operator(spin_system,{'Lz'},{electrons(2)});
    S1zS2z=operator(spin_system,{'Lz','Lz'},{electrons(1),electrons(2)});
    S1zS2p=operator(spin_system,{'Lz','L+'},{electrons(1),electrons(2)});
    S1pS2z=operator(spin_system,{'L+','Lz'},{electrons(1),electrons(2)});
    
    % Form component operators
    O11=0.5*(+S1z+S2z)+S1zS2z; O12=0.5*S2p+S1zS2p;
    O22=0.5*(+S1z-S2z)-S1zS2z; O13=0.5*S1p+S1pS2z;
    O33=0.5*(-S1z+S2z)-S1zS2z; O24=0.5*S1p-S1pS2z;
    O44=0.5*(-S1z-S2z)+S1zS2z; O34=0.5*S2p-S1zS2p;
        
    % Build the electron part of the relaxation superoperator
    R=-0.25*r1e*(O12*O12' + O12'*O12 - O11*O11 - O22*O22 + O13*O13' + O13'*O13 - O11*O11 - O33*O33 +...
                 O24*O24' + O24'*O24 - O22*O22 - O44*O44 + O34*O34' + O34'*O34 - O33*O33 - O44*O44)+...
       0.5*(r2e*(O11*O22  + O22*O11  + O11*O33 + O33*O11 + O44*O22  + O22*O44  + O33*O44 + O44*O33)+...
            r2n*(O11*O44  + O44*O11  + O22*O33 + O33*O22));
        
    % Find nuclei
    nuclei=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
        
    % Loop over nuclei
    for n=nuclei
            
        % Get basic nuclear operators
        Iz=operator(spin_system,{'Lz'},{n});
        Ip=operator(spin_system,{'L+'},{n});
        
        % Build the nuclear part of the relaxation superoperator
        R=R-0.25*r1n*(Ip*Ip'+Ip'*Ip)-r2n*Iz*Iz;
        
    end
    
    % Inform the user
    report(spin_system,'Nottingham DNP relaxation theory terms have been incorporated.');
        
end

% Add scalar relaxation of the first kind
if ismember('SRFK',spin_system.rlx.theories)

    % Inform the user
    report(spin_system,'scalar relaxation of the first kind...');
    
    % Set local assumptions
    spin_system_local=assume(spin_system,'labframe');
    
    % Get the background Hamiltonian
    H0=hamiltonian(spin_system_local);
    
    % Set local assumptions
    spin_system_local=assume(spin_system,'labframe','couplings');
    
    % Replace couplings with their modulation depths
    [rows,cols]=find(~cellfun(@isempty,spin_system.rlx.srfk_mdepth));
    spin_system_local.inter.coupling.matrix=cell(size(spin_system.inter.coupling.matrix));
    for n=1:numel(rows)
        spin_system_local.inter.coupling.matrix{rows(n),cols(n)}=...
             2*pi*eye(3)*spin_system.rlx.srfk_mdepth{rows(n),cols(n)};
    end
    
    % Get the modulated coupling Hamiltonian
    H1=hamiltonian(spin_system_local);
    
    % Erase the modifications
    clear('spin_system_local');
    
    % Eliminate numerical noise
    H0=(H0+H0')/2; H1=(H1+H1')/2;
    
    % Call the scalar relaxation superoperator
    report(spin_system,'integrating the SRFK component...');
    RS=rlx_scalar(spin_system,H0,H1,spin_system.rlx.srfk_tau_c);
    
    % Catch abuse of Redfield theory
    max_rate=max(abs(diag(RS)));
    for n=1:numel(spin_system.rlx.srfk_tau_c)
        if (1/max_rate)<spin_system.rlx.srfk_tau_c{n}(2)
            report(spin_system,['1/max(diag(R)) = ' num2str(1/max_rate) ...
                                ', srfk_tau_c = ' num2str(spin_system.rlx.srfk_tau_c{n}(2))]);
            error('T1,2>>srfk_tau_c validity condition violation in Redfield theory');
        end
    end
    
    % Add to total and tidy up 
    R=R+RS; clear('H0','H1','RS');
    
    % Inform the user
    report(spin_system,'SRFK terms have been incorporated.');
    
end

% Add scalar relaxation of the second kind
if ismember('SRSK',spin_system.rlx.theories)
    
    % Inform the user
    report(spin_system,'scalar relaxation of the second kind...');
    
    % Preallocate the answers
    R1=zeros(spin_system.comp.nspins,1); R2=zeros(spin_system.comp.nspins,1);

    % Loop over sources
    for k=spin_system.rlx.srsk_sources

        % Relaxation rates of the source spin
        Lz_k=state(spin_system,{'Lz'},{k}); Lz_k=Lz_k/norm(Lz_k,2);
        Lp_k=state(spin_system,{'L+'},{k}); Lp_k=Lp_k/norm(Lp_k,2);
        T1k=-1/real(Lz_k'*R*Lz_k); T2k=-1/real(Lp_k'*R*Lp_k);

        % Source spin quantum number
        Sk=(spin_system.comp.mults(k)-1)/2;
        
        % Inform the user
        report(spin_system, ['SRSK from spin ' num2str(k) ' (' spin_system.comp.isotopes{k} ...
                             ', T1=' num2str(T1k) ' s, T2=' num2str(T2k) ' s):']);
    
        % Loop over potential destinations
        for n=setdiff(1:spin_system.comp.nspins,...
                      spin_system.rlx.srsk_sources)

            % Only process heteronuclear pairs
            if ~strcmp(spin_system.comp.isotopes{n},spin_system.comp.isotopes{k})
            
                % Determine the coupling (rad/s)
                A=trace(spin_system.inter.coupling.matrix{n,k})/3+...
                  trace(spin_system.inter.coupling.matrix{k,n})/3;

                % Catch out-of-scope situations
                if ((2*pi/T1k)<10*abs(A))||((2*pi/T2k)<10*abs(A))
                    error('SRSK theory is not applicable: source spin relaxation is too slow');
                end
                
                % Find frequency difference
                delta_omega=spin_system.inter.basefrqs(n)-...
                            spin_system.inter.basefrqs(k);

                % Compute Abragam's expressions
                R1add=(2/3)*(A^2)*Sk*(Sk+1)*(T2k/(1+delta_omega^2*T2k^2));
                R2add=(1/3)*(A^2)*Sk*(Sk+1)*(T1k+T2k/(1+delta_omega^2*T2k^2));
                R1(n)=R1(n)+R1add; R2(n)=R2(n)+R2add;

                % Inform the user of significant terms
                if (R1add>spin_system.tols.liouv_zero)||...
                   (R2add>spin_system.tols.liouv_zero)
                    report(spin_system, ['  > onto spin ' num2str(n) ...
                                         ' (' spin_system.comp.isotopes{n}...
                                         '): +' num2str(R1add) ...
                                         ' Hz to R1, +' num2str(R2add) ' Hz to R2']);
                end
               
            end
            
        end 
    end
    
    % Clone the spin system structure
    spin_system_local=spin_system;
    spin_system_local.rlx.theories={'t1_t2'};
    spin_system_local.rlx.r1_rates=num2cell(R1');
    spin_system_local.rlx.r2_rates=num2cell(R2');
    
    % Issue a recursive call for SRSK
    report(spin_system,'recursive call for SRSK relaxation terms...'); 
    R=R+relaxation(spin_system_local);
    
    % Inform the user
    report(spin_system,'SRSK relaxation theory terms have been incorporated.');
    
end

% Print matrix density statistics
report(spin_system,['full relaxation superoperator density ' num2str(100*nnz(R)/numel(R))...
                    '%, nnz ' num2str(nnz(R)) ', sparsity ' num2str(issparse(R))]);

% Decide the fate of dynamic frequency shifts
switch spin_system.rlx.dfs
    
    case 'keep'
        
        % Do nothing and inform the user
        report(spin_system,'dynamic frequency shifts have been kept.');
        
    case 'ignore'
        
        % Kill the dynamic frequency shifts
        R=real(R);
        
        % Inform the user
        report(spin_system,'dynamic frequency shifts have been ignored.');
        
    otherwise
        
        % Complain and bomb out
        error('invalid value of the inter.rlx_dfs parameter.');
        
end

% Decide the fate of the off-diagonal components
switch spin_system.rlx.keep
    
    case 'diagonal'
        
        % Pull out the diagonal
        R=diag(diag(R));

        % Still make sure the unit state is not damped
        U=unit_state(spin_system); R=R-(U'*R*U)*(U*U');
        
        % Inform the user
        report(spin_system,'all cross-relaxation terms have been ignored.');
        
    case 'kite'
        
        % Refuse to process inappropriate cases
        if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
            error('kite option is only available for sphten-liouv formalism.');
        end
        
        % Compile the index of all longitudinal spin orders
        [~,M]=lin2lm(spin_system.bas.basis);
        long_states=find(sum(abs(M),2)==0);
        
        % Index the relaxation superoperator
        [rows,cols,vals]=find(R);
        
        % Keep self-relaxation and longitudinal cross-relaxation terms
        vals=vals.*((ismember(rows,long_states)&ismember(cols,long_states))|(rows==cols));
        
        % Recompose the relaxation superoperator
        R=sparse(rows,cols,vals,length(R),length(R)); clear('rows','cols','vals');
        
        % Inform the user
        report(spin_system,'transverse cross-relaxation terms have been ignored.');
        
    case 'secular'
        
        % Refuse to process inappropriate cases
        if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
            error('secular option is only available for sphten-liouv formalism.');
        end
        
        % Compute base frequencies of basis states
        [~,M]=lin2lm(spin_system.bas.basis);
        frequencies=sum(spin_system.inter.basefrqs.*M,2);
        
        % Index the relaxation superoperator
        [rows,cols,vals]=find(R);
        
        % Set the initial keep mask
        keep_mask=false(size(vals));
        
        % Loop over unique frequencies
        for omega=unique(frequencies)'
            
            % Find the states having the current frequency
            current_frq_group=find(frequencies==omega);
            
            % Update the keep mask
            keep_mask=keep_mask|(ismember(rows,current_frq_group)&...
                                 ismember(cols,current_frq_group));
            
        end
        
        % Recompose the relaxation superoperator
        R=sparse(rows,cols,vals.*keep_mask,length(R),length(R));
                
        % Inform the user
        report(spin_system,'secular approximation wrt Zeeman Hamiltonian...');
        report(spin_system,'non-secular cross-relaxation terms have been ignored.');
        
        % Do garbage collection
        clear('rows','cols','vals');
        
    case 'labframe'
        
        % Inform the user
        report(spin_system,'returning complete relaxation superoperator (lab frame simulations only).');
        
end

% Add non-selective damping terms
if ismember('damp',spin_system.rlx.theories)

    % Anisotropic damping option
    if exist('euler_angles','var')

        % Update the user
        report(spin_system,'adding anisotropic non-selective damping terms:');

        % Compute orientation ort (this matches alphas=0 of two-angle grids)
        ort=[0 0 1]*euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3));

        % Get the rate at the current orientation
        rate=ort*spin_system.rlx.damp_rate*ort';

    else

        % Update the user
        report(spin_system,'adding isotropic non-selective damping terms:');
        if ~isscalar(spin_system.rlx.damp_rate)
            report(spin_system,'   WARNING - damping rate anisotropy ignored');
        end
        
        % Get the isotropic part of the rate
        rate=mean(diag(spin_system.rlx.damp_rate));

    end
        
    % Build the relaxation matrix
    switch spin_system.bas.formalism
        
        case {'zeeman-hilb'}
            
            % Update the user
            report(spin_system,'   Hilbert space - unit state will be damped');
    
            % Damp everything 
            R=R-(rate/2)*unit_oper(spin_system);
            
        case {'zeeman-liouv','sphten-liouv'}
            
            % Update the user
            report(spin_system,'   Liouville space - unit state will not be damped');

            % Damp everything except unit state
            RD=-rate*unit_oper(spin_system);
            U=unit_state(spin_system);
            R=R+RD-(U'*RD*U)*(U*U');

    end
    
end
    
% Choose the thermalisation model
switch spin_system.rlx.equilibrium
    
    case 'zero'
        
        % Do nothing and inform the user
        report(spin_system,'thermalisation method: none, relaxation to unit state');
        
    case 'IME'
        
        % Inform the user
        report(spin_system,'thermalisation method: inhomogeneous master equation');
        
        % Get the equilibrium state
        if exist('euler_angles','var')
            
            % Get the Hamiltonian
            report(spin_system,'getting lab frame Hamiltonian...');
            [H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
            
            % Get the equilibrium state
            report(spin_system,'getting the equilibrium state at Hamiltonian orientation:');
            report(spin_system,['  alpha=' num2str(euler_angles(1)) ...
                                ', beta='  num2str(euler_angles(2)) ...
                                ', gamma=' num2str(euler_angles(3)) '...']);
            rho_eq=equilibrium(spin_system,H,Q,euler_angles);
            
        else
            
            % Get the Hamiltonian
            report(spin_system,'getting lab frame Hamiltonian...');
            H=hamiltonian(assume(spin_system,'labframe'),'left');
            
            % Get the equilibrium state
            report(spin_system,'getting the equilibrium state using isotropic Hamiltonian...');
            rho_eq=equilibrium(spin_system,H);
            
        end
        
        % Call the thermaliser
        R=thermalize(spin_system,R,[],[],rho_eq,'IME');
        
    case 'dibari'
        
        % Inform the user
        report(spin_system,'thermalisation method: DiBari-Levitt');
        
        % Use DiBari-Levitt thermalization
        if exist('euler_angles','var')
            
            % Get the left side Hamiltonian product superoperator
            [H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
            
            % Call the thermaliser
            R=thermalize(spin_system,R,H+orientation(Q,euler_angles),...
                         spin_system.rlx.temperature,[],'dibari');
            
            % Inform the user
            report(spin_system,'equilibrium state was computed using anisotropic Hamiltonian at');
            report(spin_system,['  alpha=' num2str(euler_angles(1)) ...
                                ', beta=' num2str(euler_angles(2)) ...
                                ', gamma=' num2str(euler_angles(3))]);
                            
        else
            
            % Get the left side Hamiltonian product superoperator
            H=hamiltonian(assume(spin_system,'labframe'),'left');
            
            % Call the thermaliser
            R=thermalize(spin_system,R,H,spin_system.rlx.temperature,[],'dibari');
            
            % Inform the user
            report(spin_system,'equilibrium state was computed using isotropic Hamiltonian.');
            
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown equilibrium specification.');
        
end

% Perform final clean-up and make array complex for later
R=complex(clean_up(spin_system,R,spin_system.tols.rlx_zero));

% Print matrix diagnostics
report(spin_system,['final relaxation superoperator density ' num2str(100*nnz(R)/numel(R))...
                    '%, nnz(R)=' num2str(nnz(R)) ', sparsity ' num2str(issparse(R))]);
R_whos=whos('R'); report(spin_system,['memory footprint of the relaxation superoperator: '...
                                      num2str(R_whos.bytes/1024^3) ' GB']);

end

% Default option values
function spin_system=defaults(spin_system)
if ~isfield(spin_system.rlx,'equilibrium')
    report(spin_system,'relaxation destination not specified, assuming zero.');
    spin_system.rlx.equilibrium='zero';
end
if ~isfield(spin_system.rlx,'dfs')
    report(spin_system,'dynamic frequency shift policy not specified, DFS will be ignored.');
    spin_system.rlx.dfs='ignore';
end
end

% Consistency enforcement
function grumble(spin_system)
if ~isfield(spin_system,'rlx')
    error('relaxation data (.rlx) is missing from the spin_system structure.');
end
if ~isfield(spin_system.rlx,'theories')
    error('relaxation data (.rlx.theories) is missing from the spin_system structure.');
end
if ~iscell(spin_system.rlx.theories)||any(~cellfun(@ischar,spin_system.rlx.theories))
    error('spin_system.rlx.theories must be a cell array of character strings.');
end
for n=1:numel(spin_system.rlx.theories)
    if ~ismember(spin_system.rlx.theories{n},{'none','damp','t1_t2','SRSK','SRFK'...
                                              'redfield','lindblad','nottingham','weizmann'})
        error('unrecognised relaxation theory specification.');
    end
end
if ( ismember('t1_t2',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv'}))
    error('extended T1,T2 relaxation theory is only available for sphten-liouv formalism.');
end
if ( ismember('redfield',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Redfield relaxation theory is only available in Liouville space.');
end
if ( ismember('lindblad',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Lindblad relaxation theory is only available in Liouville space.');
end
if ( ismember('nottingham',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Nottingham relaxation theory is only available in Liouville space.');
end
if ( ismember('weizmann',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Weizmann relaxation theory is only available in Liouville space.');
end
if ( ismember('SRFK',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('SRFK relaxation theory is only available in Liouville space.');
end
if ( ismember('SRSK',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('SRSK relaxation theory is only available in Liouville space.');
end
if ( strcmp(spin_system.rlx.equilibrium,'IME'))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('IME thermalisation is only available in Liouville space.');
end
if (strcmp(spin_system.rlx.equilibrium,'dibari'))&&...
   (spin_system.rlx.temperature==0)
    error('DiBari-Levitt thermalization cannot be used with the high-temperature approximation.')
end
end

% There are horrible people who, instead of solving a problem, tangle it up
% and make it harder to solve for anyone who wants to deal with it. Whoever
% does not know how to hit the nail on the head should be asked not to hit
% it at all.
%
% Friedrich Nietzsche

