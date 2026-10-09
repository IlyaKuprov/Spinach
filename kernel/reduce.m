% Symmetry and trajectory-level state space reduction. Tries all 
% applicable enabled reduction methods and returns a cell array
% of projectors into independently evolving reduced subspaces.
% Syntax:
%
%              projectors=reduce(spin_system,L,rho)
%
% Parameters:
%
%     L   -  Liouvillian matrix; in the compiled substance space,
%            cross-substance blocks require declared reaction records
%
%     rho -  initial state (source state screening) or
%            destination state (destination state screening)
%
% Outputs:
%
%   projectors - a cell array of projectors into independently
%                evolving reduced subspaces. The projectors are
%                to be used as follows:
%
%                    L_reduced=P'*L*P;    (for matrices)
%                    rho_reduced=P'*rho;  (for state vectors)
%
% Notes: further information on what this function does is avai-
%        lable in our papers on this subject
%
%              http://dx.doi.org/10.1016/j.jmr.2008.08.008
%              http://dx.doi.org/10.1063/1.3398146
%              http://dx.doi.org/10.1016/j.jmr.2011.03.010
%
%        Briefly, the function tries symmetry factorisation, fol-
%        lowed by zero track elimination when sys.enable contains
%        'zte', then disconnected subspace identification by path
%        tracing. Reaction records bypass spin-only symmetry blocks,
%        because atom transport need not preserve those subspaces.
%
% ilya.kuprov@weizmann.ac.il
% matthew.krzystyniak@oerc.ox.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=reduce.m>

function projectors=reduce(spin_system,L,rho)

% Check the input
grumble(spin_system,L,rho);

% Check for blanket ban
if ismember('trajlevel',spin_system.sys.disable)
    report(spin_system,'WARNING - trajectory level state space restriction disabled by user.');
    projectors{1}=1; return
end

% Reaction maps need not preserve the spin-only substance irreps
irr_projectors={}; irr_dimensions=[];
if size(L,1)==spin_system.bas.offsets(end)&&isempty(spin_system.chem.reactions)
    for s=1:spin_system.bas.nsubst
        for n=1:numel(spin_system.bas.sym_fact(s).irr_dimensions)
            P=spin_system.bas.sym_fact(s).irr_projectors{n};
            irr_projectors{end+1}=[sparse(spin_system.bas.offsets(s),size(P,2)); P;...
                sparse(spin_system.bas.offsets(end)-spin_system.bas.offsets(s+1),size(P,2))]; %#ok<AGROW>
        end
        irr_dimensions=[irr_dimensions; spin_system.bas.sym_fact(s).irr_dimensions(:)]; %#ok<AGROW>
    end
end

% Identify unit directions only in the compiled spherical-tensor space
unit_states=sparse(size(L,1),0);
if strcmp(spin_system.bas.formalism,'sphten-liouv')&&...
   (size(L,1)==spin_system.bas.offsets(end))
    unit_states=sparse(spin_system.bas.offsets(1:end-1)+1,...
                       1:spin_system.bas.nsubst,1,size(L,1),spin_system.bas.nsubst);
end

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'zeeman-hilb'
        
        % If a cell array is supplied, make a representative density matrix
        if iscell(rho)
            rho_rep=abs(rho{1});
            for n=2:numel(rho)
                rho_rep=rho_rep+abs(rho{n});
            end
            rho=rho_rep/numel(rho);
        end
        
        % Run symmetry factorization
        if ismember('symmetry',spin_system.sys.disable)
            
            % Issue a reminder to the user
            report(spin_system,'WARNING - permutation symmetry treatment disabled by user.');
            
            % Return unit matrix
            projectors{1}=1;
            
        elseif isempty(irr_projectors)
            
            % Inform the user
            report(spin_system,'no permutation symmetry information has been supplied.');
            
            % Return unit matrix
            projectors{1}=1;
            
        else
            
            % Decide which irreps to keep
            n_irreps=numel(irr_projectors);
            irrep_keep_index=true(n_irreps,1);
            for n=1:n_irreps
                
                % Check the irrep contribution to the total norm
                if irr_dimensions(n)==0
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ...
                                        ' has no states in it - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                elseif norm(irr_projectors{n}'*rho*... %#NORMOK
                            irr_projectors{n},1)<spin_system.tols.irrep_drop
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ' has less than '...
                                        num2str(spin_system.tols.irrep_drop) ...
                                        ' of the initial state norm - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                else
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ' contains '...
                                        num2str(irr_dimensions(n)) ...
                                        ' states - kept.']);
                    
                end
                
            end
            
            % Compile the projector array
            projectors=irr_projectors(irrep_keep_index);
        
        end
        
    case 'zeeman-wavef'

        % If a stack is supplied, choose a representative wavefunction
        if size(rho,2)>1, rho=mean(abs(rho),2); end

        % Run symmetry factorization
        if ismember('symmetry',spin_system.sys.disable)

            % Issue a reminder to the user
            report(spin_system,'WARNING - permutation symmetry treatment disabled by user.');

            % Return unit matrix
            projectors{1}=1;

        elseif isempty(irr_projectors)

            % Inform the user
            report(spin_system,'no permutation symmetry information has been supplied.');

            % Return unit matrix
            projectors{1}=1;

        else

            % Decide which irreps to keep
            n_irreps=numel(irr_projectors);
            irrep_keep_index=true(n_irreps,1);
            for n=1:n_irreps

                % Check the irrep contribution to the total norm
                if irr_dimensions(n)==0

                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ' has dimension zero - dropped.']);

                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;

                elseif norm(irr_projectors{n}'*rho,1)<spin_system.tols.irrep_drop

                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(irr_dimensions(n)) ' has less than '...
                                        num2str(spin_system.tols.irrep_drop) ' of the state norm - dropped.']);

                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;

                else

                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(irr_dimensions(n)) ' is active - kept.']);

                end

            end

            % Compile the projector array
            projectors=irr_projectors(irrep_keep_index);

        end

    case {'zeeman-liouv','sphten-liouv'}

        % If a stack is supplied, choose a representative state vector
        if size(rho,2)>1, rho=mean(abs(rho),2); end
        
        % Run symmetry factorization and ZTE
        if ismember('symmetry',spin_system.sys.disable)
            
            % Issue a reminder to the user
            report(spin_system,'WARNING - permutation symmetry treatment disabled by user.');
            
            % Run zero track elimination
            report(spin_system,'attempting zero track elimination...');
            projectors{1}=zte(spin_system,L,rho);
            
        elseif isempty(irr_projectors)
            
            % Inform the user
            report(spin_system,'no independent permutation symmetry blocks are available.');
            
            % Run zero track elimination
            report(spin_system,'attempting zero track elimination...');
            projectors{1}=zte(spin_system,L,rho);
            
        else
            
            % Decide which irreps to keep
            n_irreps=numel(irr_projectors); irrep_keep_index=true(n_irreps,1);
            for n=1:n_irreps
                
                % Check the irrep contribution to the total norm
                if irr_dimensions(n)==0
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ' has dimension zero - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                elseif norm(irr_projectors{n}'*rho,1)<spin_system.tols.irrep_drop&&...
                       (nnz(irr_projectors{n}'*unit_states)==0)
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(irr_dimensions(n)) ' has less than '...
                                        num2str(spin_system.tols.irrep_drop) ' of the state norm - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                else
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(irr_dimensions(n)) ' is active - kept.']);
                    
                end
                
            end
            
            % Compile the projector array
            projectors=irr_projectors(irrep_keep_index);
            
            % Loop over permutation group irreps
            for n=1:numel(projectors)
                
                % Report to the user
                report(spin_system,['irrep #' num2str(n) ', attempting zero track elimination...']);
                
                % Mark unit support in the projected coordinates, not the original offsets
                local_system=spin_system;
                unit_rows=find(any(projectors{n}'*unit_states,2));
                local_system.bas.offsets=[unit_rows-1; size(projectors{n},2)];

                % Run zero track elimination
                zte_projector=zte(local_system,projectors{n}'*L*projectors{n},projectors{n}'*rho);
                
                % Project the projectors
                projectors{n}=projectors{n}*zte_projector;
                
            end
            
        end
        
        % Run path tracing
        for n=1:numel(projectors)
            
            % Inform the user
            report(spin_system,['path-tracing subspace #' num2str(n) '...']);
            
            % Include unit support in population screening after projection
            rho_screen=abs(projectors{n}'*rho)+sum(abs(projectors{n}'*unit_states),2);

            % Run the path tracing
            pt_projectors=path_trace(spin_system,projectors{n}'*L*projectors{n},rho_screen);
            
            % Project the projectors
            for k=1:numel(pt_projectors)
                pt_projectors{k}=projectors{n}*pt_projectors{k};
            end
            projectors{n}=pt_projectors(:); %#ok<AGROW>
            
        end
        
        % Flatten out the cell array
        projectors=vertcat(projectors{:});
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% Consistency enforcement
function grumble(spin_system,L,rho) %#ok<INUSD,INUSL>
if ~isnumeric(L)
    error('L must be numeric.');
end
if size(L,1)~=size(L,2)
    error('L must be a square matrix.');
end
if spin_system.bas.nsubst>1&&size(L,1)==spin_system.bas.offsets(end)&&...
   isempty(spin_system.chem.reactions)
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        if nnz(L(idx,:))~=nnz(L(idx,idx))
            error('Spinach:reduce:crossSubstanceGenerator',...
                  'L must not contain cross-substance blocks (substance %d).',n);
        end
    end
end
end

% For centuries, the battle of morality was fought between those who claimed
% that your life belongs to God and those who claimed that it belongs to your
% neighbours - between those who preached that the good is self-sacrifice for
% the sake of ghosts in heaven and those who preached that the good is self-
% sacrifice for the sake of incompetents on earth. And no one came to say that
% your life belongs to you and the good is to live it.
%
% Ayn Rand, "Atlas Shrugged"

