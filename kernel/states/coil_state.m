% Generates Hilbert space density matrices and Liouville space state 
% vectors from their human-readable descriptions. Syntax:
%
%              rho=coil_state(spin_system,states,spins,method)
%
% Parameters:
%
% 1. If states is a string and spins is a string
%
%                      states='Lz'; spins='13C';
%
% the function returns the sum of the corresponding single-spin densi-
% ty matrices (Hilbert space) or state vectors (Liouville space) on 
% all spins of that type. Valid labels for states in this type of call
% are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreduci-
% ble spherical tensor, l and m are integers), 'CTx', 'CTy', 'CTz',
% 'CT+','CT-' (central transition operators in the Zeeman basis). Va-
% lid labels for spins are standard isotope names, as well as 'elect-
% rons', 'nuclei', and 'all'.
%
% 2. If states is a string and spins is a vector
%
%                      states='Lz'; spins=[1 2 4];
%
% the function returns the sum of all single-spin density matrices
% (Hilbert space) or state vectors (Liouville space) for all spins
% with the specified numbers. Valid labels for states are the same as
% in Item 1 above.
%
% 3. If states is a cell array of strings and spins is a cell array
%    of numbers:
%
%                    states={'Lz','L+'}; spins={1,2};
%
% then a product state density matrix (Hilbert space) or state vector
% (Liouville space) is produced. In the case above, Spinach will gene-
% rate LzS+ density matrix in Hilbert space or its state vector in Li-
% ouville space. Valid labels for operators are the same as in Item 1
% above.
%
% 4. For wavefunction formalism, states must be specified as an array
%    of projection quantum numbers on all spins; pass an empty spin list,
%    for example, in a {'1H','1H','14N'} system:
%
%                   psi=coil_state(spin_system,[-1/2 1/2 0],[],'exact')
%
% Method argument has the following effect in sphten-liouv formalism:
%
%    'cheap'  - the state vector is generated without
%               normalisation. For very large spin sys-
%               tems this is much faster
%
%    'exact'  - exact state vector with correct normalisation
%
% All four arguments are required; method is ignored in Zeeman formalisms.
%
% The result is unweighted in every formalism. Use state for initial
% populations weighted by substance concentrations. Wavefunction output is
% a direct sum of unweighted product kets, one per substance; this storage
% representation does not support concentration weighting or chemistry.
%
% Outputs:
%
%     rho     - a Hilbert space density matrix or a Liouville
%               space state vector
%
% d.savostyanov@soton.ac.uk
% luke.edwards@ucl.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=coil_state.m>

function rho=coil_state(spin_system,states,spins,method)

% Check consistency
grumble(spin_system,states,spins,method);

% Reject product states spanning substances
if iscell(spins), which_subst(spin_system,cell2mat(spins)); end

% Preserve substance membership before identity factors lose their spin labels
if strcmp(spin_system.bas.formalism,'sphten-liouv')&&spin_system.bas.nsubst>1
    if ischar(spins)
        switch spins
            case 'all'
                spins=1:spin_system.comp.nspins;
            case 'electrons'
                spins=find(cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
            case 'nuclei'
                spins=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
            otherwise
                spins=find(strcmp(spins,spin_system.comp.isotopes));
        end
        if isempty(spins), error('no such spins in the system.'); end
    end
    if isnumeric(spins)
        rho=sparse(spin_system.bas.offsets(end),1);
        for n=1:numel(spins)
            rho=rho+coil_state(spin_system,{states},{spins(n)},method);
        end
        return;
    end
    subst=which_subst(spin_system,cell2mat(spins));
end

% Get the unit state
switch spin_system.bas.formalism

    case 'sphten-liouv'

        % Build unweighted substance unit coordinates
        unit=sparse(spin_system.bas.offsets(1:end-1)+1,1,1,...
                    spin_system.bas.offsets(end),1);
        if spin_system.bas.nsubst>1
            unit=sparse(spin_system.bas.offsets(subst)+1,1,1,...
                        spin_system.bas.offsets(end),1);
        end
        
    case 'zeeman-liouv'

        % Vectorise the unweighted Hilbert identity
        blocks=cell(spin_system.bas.nsubst,1);
        for n=1:spin_system.bas.nsubst
            block=speye(prod(spin_system.comp.mults(spin_system.chem.parts{n})));
            blocks{n}=block(:);
        end
        unit=vertcat(blocks{:});

end

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Choose the state vector generation method
        switch method
            
            % Careless normalisation
            case 'cheap'
                
                % Parse the specification
                [opspecs,coeffs]=human2opspec(spin_system,states,spins);
                
                % Allocate the unweighted direct-sum state
                rho=sparse(spin_system.bas.offsets(end),1);

                % Locate each operator in its hosting descriptor
                for n=1:numel(opspecs)
                    active_spins=find(opspecs{n});
                    if isempty(active_spins)
                        rho=rho+coeffs(n)*unit;
                        continue;
                    end
                    subst=which_subst(spin_system,active_spins);
                    descriptor=opspecs{n}(spin_system.chem.parts{subst});
                    [present,index]=ismember(descriptor,spin_system.bas.basis{subst},'rows');
                    if ~present
                        error('the requested state is not present in the basis.');
                    end
                    index=index+spin_system.bas.offsets(subst);
                    rho(index)=rho(index)+coeffs(n);
                end

            % Careful normalisation
            case 'exact'
                
                % Apply a left side product superoperator to the unit state
                rho=operator(spin_system,states,spins,'left')*unit;
            
            otherwise
                
                % Complain and bomb out
                error('unknown state generation method.');
                
        end
        
    case 'zeeman-liouv'

        % Apply a left side product superoperator to the unit state
        rho=operator(spin_system,states,spins,'left')*unit; 
                
    case 'zeeman-hilb'
        
        % Generate a Hilbert space operator
        rho=operator(spin_system,states,spins);

    case 'zeeman-wavef'

        % Store one unweighted product ket per substance
        blocks=cell(spin_system.bas.nsubst,1);
        for s=1:spin_system.bas.nsubst
            psi=1;
            for n=spin_system.chem.parts{s}

                % Select the requested local projection eigenstate
                current_mult=spin_system.comp.mults(n);
                current_spin=(current_mult-1)/2;
                levels=fliplr((-current_spin):(current_spin));
                current_psi=double(levels==states(n));
                psi=kron(psi,transpose(current_psi));

            end
            blocks{s}=psi;
        end
        rho=vertcat(blocks{:});

    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
    
end

end

% Input validation function
function grumble(spin_system,states,spins,method)

if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))
    error('basis set information is missing, run basis() before calling this function.');
end
if ~ischar(spin_system.bas.formalism)
    error('formalism specification must be a character string.');
end
if ~ismember(spin_system.bas.formalism,{'zeeman-hilb', 'zeeman-liouv',...
                                        'sphten-liouv','zeeman-wavef'})
    error('unknown formalism specification.');
end

if ~ischar(method)
    error('method must be a character string.')
elseif ~ismember(method, {'cheap', 'exact'})
    error('unknown method specification.');
end

if (~(ischar(states)&&ischar(spins)))&&...
   (~(iscell(states)&&iscell(spins)))&&...
   (~(ischar(states)&&isnumeric(spins)))&&...
   (~(isnumeric(states)&&isempty(spins)))
    error('invalid state specification.');
end
if isnumeric(states)&&(numel(states)~=spin_system.comp.nspins)
    error('numel(states) must match number of spins in the system.')
end
if isnumeric(states)&&any(mod(states,0.5)~=0,'all')
    error('spin projection numbers must be integer or half-integer.');
end
if isnumeric(states)
    spin_qn=(spin_system.comp.mults(:).'-1)/2;
    if any(abs(states(:).')>spin_qn)||any(mod(states(:).'+spin_qn,1)~=0)
        error('each projection quantum number must be an allowed level of its spin.');
    end
end
if iscell(states)&&iscell(spins)&&(numel(states)~=numel(spins))
    error('spins and operators cell arrays should have the same number of elements.');
end
if iscell(states)&&any(~cellfun(@ischar,states))
    error('all elements of the operators cell array should be strings.');
end
if isnumeric(spins)&&(~isempty(spins))
    if (~isreal(spins))||(~isrow(spins))||any(mod(spins,1)~=0)||any(spins<1)
        error('when numeric, spins must be a row of positive integers.');
    end
    if numel(spins)~=numel(unique(spins))
        error('spin list must not have any repetitions.');
    end
end
if iscell(spins)
    if isempty(spins)
        error('when a cell array, spin list cannot be empty.');
    end  
    for n=1:numel(spins)
        if (~isreal(spins{n}))||(mod(spins{n},1)~=0)||(spins{n}<1)
            error('when a cell array, spins must contain positive integers.');
        end
    end
    spins=cell2mat(spins(:));
    if numel(spins)~=numel(unique(spins))
        error('spin list must not have any repetitions.');
    end
end
end

% Aggressive public displays of virtue are where 
% the morally deplorable hide.
%
% Milo Yiannopoulos


