% Chemical kinetics superoperator. Syntax:
%
%                      K=kinetics(spin_system)
%
% Parameters:
%
%    spin_system  -  Spinach spin system description object
%                    produced as described in the spin system
%                    and basis specification sections, of the
%                    of the online manual. All adjustable pa-
%                    rameters are described in the chemical 
%                    kinetics parameters section.
%
% Outputs:
%
%    K            -  kinetics superoperator
%
% Note: a large variety of chemical reaction models is supported,
%       see the chemical kinetics parameters section of the onli-
%       ne manual.
%
% Note: Spinach context functions include kinetics superoperators
%       into the total Liovillian automatically.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% hannah.hogben@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kinetics.m>

function K=kinetics(spin_system)

% Preallocate the answer
K=mprealloc(spin_system,0);

% Find chemical reactions
[sources,destins,rates]=find(spin_system.chem.rates);

% Loop over chemical reactions
for n=1:numel(rates)
    
    % Catch incorrect calls
    if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
        error('chemical reaction treatment is only available for sphten-liouv formalism.');
    end
    
    % Get the spins involved
    source_spins=spin_system.chem.parts{sources(n)};
    destin_spins=spin_system.chem.parts{destins(n)};
    
    % Get the states involved
    source_states=logical(sum(spin_system.bas.basis(:,source_spins),2));
    destin_states=logical(sum(spin_system.bas.basis(:,destin_spins),2));
    
    % Make sure basis tables match
    if ~isequal(spin_system.bas.basis(source_states,source_spins),...
                spin_system.bas.basis(destin_states,destin_spins))
        error('spin systems on either side of the reaction arrow have different topologies or basis sets.');
    end

    % Move to integer indexing
    source_states=find(source_states);
    destin_states=find(destin_states);
    
    % Update the kinetics superoperator
    K=K+rates(n)*sparse(source_states,destin_states,ones(size(source_states)),size(K,1),size(K,2));

end

% Find magnetization fluxes
[source_spins,destin_spins,flux_rate]=find(spin_system.chem.flux_rate);

% Process magnetization fluxes
if ~isempty(flux_rate)
    
    % Catch incorrect calls
    if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
        error('magnetization flux treatment is only available for sphten-liouv formalism.');
    end
    
    % Index single- and multi-spin orders (sso and mso)
    sso_state_mask=(sum(logical(spin_system.bas.basis),2)==1);
    mso_state_mask=(sum(logical(spin_system.bas.basis),2)>1 );
    
    % Decide how to proceed
    switch spin_system.chem.flux_type
        case 'intramolecular'
            report(spin_system,'intramolecular magnetization flux, correlations will be kept.');
        case 'intermolecular'
            report(spin_system,'intermolecular magnetization flux, correlations will be damped.');
        otherwise
            error('unrecognized magnetization flux type.')
    end
    
    % Loop over the fluxes
    for n=1:numel(flux_rate)
        
        % Find single-spin sources and destinations
        sso_source_state_mask=sso_state_mask&(spin_system.bas.basis(:,source_spins(n))~=0);
        sso_destin_state_mask=sso_state_mask&(spin_system.bas.basis(:,destin_spins(n))~=0);
        
        % Identify stationary states
        sso_static_state_mask=sso_source_state_mask&sso_destin_state_mask;
        sso_source_state_mask=xor(sso_source_state_mask,sso_static_state_mask);
        sso_destin_state_mask=xor(sso_destin_state_mask,sso_static_state_mask);
        
        % Make sure subspaces match
        if ~isequal(spin_system.bas.basis(sso_source_state_mask,source_spins(n)),...
                    spin_system.bas.basis(sso_destin_state_mask,destin_spins(n)))
            error('spin systems on either side of the reaction arrow have different topology or basis sets.');
        end
        
        % Move to integer indexing
        sso_source_state_index=find(sso_source_state_mask);
        sso_destin_state_index=find(sso_destin_state_mask);
        
        % Update the kinetics superoperator
        K=K+flux_rate(n)*sparse(sso_destin_state_index,...
                                sso_source_state_index,...
                                ones(size(sso_source_state_index)),size(K,1),size(K,2));
        K=K-flux_rate(n)*sparse(sso_source_state_index,...
                                sso_source_state_index,...
                                ones(size(sso_source_state_index)),size(K,1),size(K,2));
        
        % Decide the fate of multi-spin orders
        switch spin_system.chem.flux_type
            
            case 'intramolecular' 
        
                % Find multi-spin sources and destinations
                mso_source_state_mask=mso_state_mask&(spin_system.bas.basis(:,source_spins(n))~=0);
                mso_destin_state_mask=mso_state_mask&(spin_system.bas.basis(:,destin_spins(n))~=0);
                
                % Identify stationary states
                mso_static_state_mask=mso_source_state_mask&mso_destin_state_mask;
                mso_source_state_mask=xor(mso_source_state_mask,mso_static_state_mask);
                mso_destin_state_mask=xor(mso_destin_state_mask,mso_static_state_mask);
                
                % Make sure subspaces match
                if ~isequal(spin_system.bas.basis(mso_source_state_mask,source_spins(n)),...
                            spin_system.bas.basis(mso_destin_state_mask,destin_spins(n)))
                    error('spin systems on either side of the reaction arrow have different topology or basis sets.');
                end
                
                % Move to integer indexing
                mso_source_state_index=find(mso_source_state_mask);
                mso_destin_state_index=find(mso_destin_state_mask);
                
                % Update the kinetics superoperator
                K=K+flux_rate(n)*sparse(mso_destin_state_index,...
                                        mso_source_state_index,...
                                        ones(size(mso_source_state_index)),size(K,1),size(K,2));
                K=K-flux_rate(n)*sparse(mso_source_state_index,...
                                        mso_source_state_index,...
                                        ones(size(mso_source_state_index)),size(K,1),size(K,2));
        
            case 'intermolecular'
                
                % Find multi-spin sources 
                mso_source_state_mask=mso_state_mask&(spin_system.bas.basis(:,source_spins(n))~=0);
                
                % Generate indices
                mso_source_state_index=find(mso_source_state_mask);
        
                % Update the kinetics superoperator
                K=K-flux_rate(n)*sparse(mso_source_state_index,...
                                        mso_source_state_index,...
                                        ones(size(mso_source_state_index)),size(K,1),size(K,2));
        
            otherwise
                
                % Complain and bomb out
                error('unrecognized type of magnetization flux.');
                
        end
        
    end

end

% Process radical pair recombination
if (~isempty(spin_system.chem.rp_theory))&&...
   (~strcmp(spin_system.chem.rp_theory,'off'))
    
    % Catch incorrect calls
    if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
        error('magnetization flux treatment is only available for sphten-liouv formalism.');
    end

    % Report to the user
    report(spin_system,['generating kinetics superoperator using the ' spin_system.chem.rp_theory ' model...'])
    
    % Get the singlet and triplet projector product superoperators
    unit=unit_oper(spin_system);
    LzSz_l=operator(spin_system,{'Lz','Lz'},num2cell(spin_system.chem.rp_electrons),'left');
    LzSz_r=operator(spin_system,{'Lz','Lz'},num2cell(spin_system.chem.rp_electrons),'right');
    LpSm_l=operator(spin_system,{'L+','L-'},num2cell(spin_system.chem.rp_electrons),'left');
    LpSm_r=operator(spin_system,{'L+','L-'},num2cell(spin_system.chem.rp_electrons),'right');
    LmSp_l=operator(spin_system,{'L-','L+'},num2cell(spin_system.chem.rp_electrons),'left');
    LmSp_r=operator(spin_system,{'L-','L+'},num2cell(spin_system.chem.rp_electrons),'right');
    singlet_l=unit/4-(LzSz_l+0.5*(LpSm_l+LmSp_l)); triplet_l=unit-singlet_l;
    singlet_r=unit/4-(LzSz_r+0.5*(LpSm_r+LmSp_r)); triplet_r=unit-singlet_r;
    
    % Assemble the recombination kinetics superoperator
    switch spin_system.chem.rp_theory
        
        case 'exponential'
            
            K=K-sum(spin_system.chem.rp_rates)*speye(size(K));
        
        case 'haberkorn'
            
            K=K-0.5*(spin_system.chem.rp_rates(1)*(singlet_l+singlet_r)+spin_system.chem.rp_rates(2)*(triplet_l+triplet_r));
            
        case 'jones-hore'
            
            K=K-sum(spin_system.chem.rp_rates)*unit+spin_system.chem.rp_rates(1)*triplet_l*triplet_r+spin_system.chem.rp_rates(2)*singlet_l*singlet_r;
            
        otherwise
            
            error('unknown radical pair kinetics model.');
            
    end
    
end

% Warn if no kinetics
if nnz(K)==0
    report(spin_system,'no significant kinetics specified.');
    report(spin_system,'kinetics superoperator set to zero.');
end

end

% The egoist is the absolute sense is not the man who sacrifices others. He
% is the man who stands above the need of using others in any manner. He do-
% es not function through them. He is not concerned with them in any primary
% matter. Not in his aim, not in his motive, not in his thinking, not in his
% desires, not in the source of his energy. He does not exist for any other 
% man - and he asks no other man to exist for him. This is the only form of
% brotherhood and mutual respect possible between men.
%
% Ayn Rand, "The Fountainhead"

