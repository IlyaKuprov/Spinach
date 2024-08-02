% Generates Hilbert space density matrices and Liouville space state 
% vectors from their human-readable descriptions. Syntax:
%
%              rho=state(spin_system,states,spins,method)
%
% This function supports three types of calls:
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
%                      states={'Lz','L+'}; spins={1,2};
%
% then a product state density matrix (Hilbert space) or state vector
% (Liouville space) is produced. In the case above, Spinach will gene-
% rate LzS+ density matrix in Hilbert space or its state vector in Li-
% ouville space. Valid labels for operators are the same as in Item 1
% above.
%
% Method argument has the following effect in sphten-liouv formalism:
%
%    'cheap'  - the state vector is generated without
%               normalisation. For very large spin sys-
%               tems this is much faster
%
%    'exact'  - exact state vector with correct normalisation,
%               this is the default when the last argument is
%               skipped in the function call
%
%    'chem'   - the exact state vector weighted with the 
%               concentrations specified in inter.chem.concs
%               field under chemical kinetics parameters
%
% This option is ignored in zeeman-hilb and zeeman-liouv formalisms
% because there are no cheap shortcuts and kinetics is not available.
%
% Outputs:
%
%     rho     - a Hilbert space density matrix or a Liouville
%               space state vector
% 
% d.sayostyanov@soton.ac.uk
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=state.m>

function rho=state(spin_system,states,spins,method)

% Validate the input
grumble(spin_system,states,spins);

% Default is to use consistent state norms
if ~exist('method','var'), method='exact'; end

% Get the unit state
switch spin_system.bas.formalism

    case 'sphten-liouv'

        % Unit population of T(0,0) state, normalisation is
        % such because prod(spin_system.comp.mults) can be-
        % come too large for double precision arithmetic
        unit=sparse(1,1,1,size(spin_system.bas.basis,1),1);
        
    case 'zeeman-liouv'

        % Stretched unit matrix, normalisation matched to 
        % the Hilbert space because systems are small
        unit=speye(prod(spin_system.comp.mults)); unit=unit(:);

end

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Choose the state vector generation methos
        switch method
            
            % Careless normalisation
            case 'cheap'
                
                % Parse the specification
                [opspecs,coeffs]=human2opspec(spin_system,states,spins);
                
                % Compute correlation orders
                correlation_orders=sum(logical(spin_system.bas.basis),2);
                
                % Locate each operator in the basis
                indices=zeros(size(coeffs));
                parfor n=1:numel(opspecs) %#ok<*PFBNS>
                    
                    % Find states with the same correlation order
                    possibilities=(correlation_orders==nnz(opspecs{n}));
                    
                    % Pin down the required state
                    for k=find(opspecs{n})
                        possibilities=and(possibilities,spin_system.bas.basis(:,k)==opspecs{n}(k)); 
                    end

                    % Double-check
                    if nnz(possibilities)>1
                        error('basis descriptor ambiguity detected.');
                    elseif nnz(possibilities)<1
                        error('the requested state is not present in the basis.');
                    end
                    
                    % Locate the state 
                    indices(n)=find(possibilities);
                    
                end
                
                % Assemble the state vector
                nrows=size(spin_system.bas.basis,1); ncols=1;
                rho=sparse(indices,ones(size(indices)),coeffs,nrows,ncols);

            % Careful normalisation
            case 'exact'
                
                % Apply a left side product superoperator to the unit state
                rho=operator(spin_system,states,spins,'left')*unit;
            
            % Chemical weighing
            case 'chem'
                
                % Parse the specification
                [opspecs,coeffs]=human2opspec(spin_system,states,spins);
                
                % Preallocate the state vector
                rho=spalloc(size(spin_system.bas.basis,1),1,0);

                % Get the basis dimension
                matrix_dim=size(spin_system.bas.basis,1);
                
                % Sum the states with concentrations
                for n=1:numel(opspecs)
                    
                    % Identify active spins
                    active_spins=find(opspecs{n});
                    
                    % Find out which chemical species they are in
                    species=true(1,numel(spin_system.chem.parts)); 
                    for k=1:numel(active_spins)
                        species=species&cellfun(@(x)ismember(active_spins(k),x),spin_system.chem.parts);
                    end
                    
                    % Check state validity
                    if nnz(species)~=1
                        error('the spin state requested crosses chemical species boundaries.');
                    end
                    
                    % Adjust the coefficient
                    coeffs(n)=coeffs(n)*spin_system.chem.concs(species);
                    
                    % Get the operator
                    A=p_superop(spin_system,opspecs{n},'left');
                    A=sparse(A(:,1),A(:,2),A(:,3),matrix_dim,matrix_dim);

                    % Get the state vector
                    rho=rho+coeffs(n)*A*unit;
                    
                end
                
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
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
    
end

end

% Input validation function
function grumble(spin_system,states,spins)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~(ischar(states)&&ischar(spins)))&&...
   (~(iscell(states)&&iscell(spins)))&&...
   (~(ischar(states)&&isnumeric(spins)))
    error('invalid state specification.');
end
if iscell(states)&&iscell(spins)&&(numel(states)~=numel(spins))
    error('spins and operators cell arrays should have the same number of elements.');
end
if iscell(states)&&any(~cellfun(@ischar,states))
    error('all elements of the operators cell array should be strings.');
end
if isnumeric(spins)
    if (~isreal(spins))||(~isrow(spins))||any(mod(spins,1)~=0)||any(spins<1)
        error('when numeric, spins must be a row of positive integers.');
    end
    if numel(spins)~=numel(unique(spins))
        error('spin list must not have any repetitions.');
    end
end
if iscell(spins)
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
if isempty(spins), error('spin list cannot be empty.'); end
end

% Aggressive public displays of virtue are where 
% the morally deplorable hide.
%
% Milo Yiannopoulos

