% Generates Hilbert space operators or Liouville space superoperators from
% their human-readable descriptions. Syntax:
%
%        A=operator(spin_system,operators,spins,operator_type,format)
%
% This function supports three types of calls:
%
% 1. If operators is a string and spins is a string
%
%	                    operators='Lz'; spins='13C';
%
% the function returns the sum of the corresponding single-spin operators 
% (Hilbert space) or superoperators (Liouville space) on all spins of that
% type. Valid labels for states in this type of call are 'E' (identity),
% 'Lz', 'Lx', 'Ly', 'L+', 'L-', 'Tl,m' (irreducible spherical tensor, l 
% and m are integers), 'CTx', 'CTy', 'CTz', 'CT+', 'CT-' (central transi-
% tion operators in the Zeeman basis). Valid labels for spins are standard 
% isotope names, as well as 'electrons', 'nuclei', and 'all'.
%
% 2. If operators is a string and spins is a vector
%
%                      operators='Lz'; spins=[1 2 4];
%
% the function returns the sum of all single-spin  operators (Hilbert space)
% or superoperators (Liouville space) for all spins with the specified num-
% bers. Valid labels for operators are the same as in Item 1 above.
%
% 3. If operators is a cell array of strings and spins is a cell array of
% numbers
%
%                     operators={'Lz','L+'}; spins={1,2};
%
% then a product operator (Hilbert space) or its superoperator (Liouville
% space) is produced. In the case above, Spinach will generate LzS+ in Hil-
% bert space or its specified superoperator in Liouville space. Valid la-
% bels for operators are the same as in Item 1 above.
%
% In Liouville space calculations, operator_type can be set to:
%
%            'left' - produces left side product superoperator
%
%           'right' - produces right side product superoperator
%
%            'comm' - produces commutation superoperator (default)
%
%           'acomm' - produces anticommutation superoperator
%
% In Hilbert space calculations operator_type parameter is ignored, and the
% operator itself is always returned.
%
% The format parameter refers to the format of the output: 'csc' returns a
% Matlab sparse matrix, 'xyz' returns a [rows, cols, vals] array.
%
% Outputs:
%
%    A   - a CSC sparse (default) or a [rows, cols, vals] repre-
%          sentation of a spin operator or superoperator.
%
% Notes: WARNING - a product of two commutation superoperators is NOT a com-
%        mutation superoperator of a product. In Liouville space, you cannot
%        generate single-spin superoperators and multiply them up.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=operator.m>

function A=operator(spin_system,operators,spins,operator_type,format)

% The default superoperator type is commutation
if ~exist('operator_type','var'), operator_type='comm'; end

% The default sparse matrix format is CSC
if ~exist('format','var'), format='csc'; end

% Check consistency
grumble(spin_system,operators,spins,operator_type,format); tic;

% Check disk cache
if ismember('op_cache',spin_system.sys.enable)

    % Combine specification, isotopes, and basis hash
    op_hash=md5_hash({operators,spins,operator_type,...
                      format,spin_system.comp.iso_hash,...
                      spin_system.bas.basis_hash});

    % Generate the cache record name in the global scratch (for later reuse)
    filename=[spin_system.sys.scratch filesep 'spinach_op_' op_hash '.mat'];

    % Load the operator from the cache record
    if exist(filename,'file'), load(filename,'A'); return; end
    
end

% Parse the human specification into Spinach notation
[opspecs,coeffs]=human2opspec(spin_system,operators,spins);

% Start a cell array
A=cell(numel(opspecs),1);

% Decide how to proceed
switch spin_system.bas.formalism
    
    % Spherical tensor basis
    case 'sphten-liouv'
        
        % Build summation terms
        parfor n=1:numel(opspecs)

            % Check physics type
            if all(opspecs{n}>=0,'all')

                % For spins: get the superoperator
                A{n}=p_superop(spin_system,opspecs{n},operator_type);

                % For spins: apply the coefficient
                A{n}(:,3)=coeffs(n)*A{n}(:,3);

            else

                % For cavities and phonons: complain and bomb out
                error('Cavities and phonons not yet available in sphten-liouv.');

            end

        end
                
    % Zeeman basis formalisms
    case {'zeeman-hilb','zeeman-liouv'}

        % Parallelisation efficiency
        mults=spin_system.comp.mults;
        formalism=spin_system.bas.formalism;

        % Build summation terms
        parfor n=1:numel(opspecs)

            % Start the kron
            B=sparse(coeffs(n));
            
            % Over particles
            for k=1:numel(opspecs{n})

                % Check physics type
                if opspecs{n}(k)>=0

                    % Spin: get spin state indices
                    [L,M]=lin2lm(opspecs{n}(k));

                    % Spin: get irreducible spherical tensors
                    ist=irr_sph_ten(mults(k),L); %#ok<PFBNS>

                    % Spin: get operator
                    T=ist{L-M+1};

                else

                    % Cavities and phonons: get operator
                    T=boson_oper(mults(k),opspecs{n}(k));

                end
                
                % Update the kron
                B=kron(B,T);

            end
            
            % Move to Liouville space if necessary
            if strcmp(formalism,'zeeman-liouv')
                B=hilb2liouv(B,operator_type);
            end

            % Convert sparse array from CSC to XYZ indexing
            [rows,cols,vals]=find(B); A{n}=[rows cols vals];

        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism.');
        
end

% XYZ format sum
A=cell2mat(A);

% Convert to CSC format
if strcmp(format,'csc')

    % Decide operator dimension
    switch spin_system.bas.formalism

        case 'sphten-liouv'

            % As per the basis set specification
            matrix_dim=size(spin_system.bas.basis,1);

        case 'zeeman-hilb'
            
            % Entire Hilbert space
            matrix_dim=prod(spin_system.comp.mults);

        case 'zeeman-liouv'
            
            % Entire Liouville space
            matrix_dim=prod(spin_system.comp.mults)^2;

        otherwise
        
            % Complain and bomb out
            error('unknown formalism.');

    end

    % Make a sparse matrix, making sure it's complex for later
    A=sparse(A(:,1),A(:,2),complex(A(:,3)),matrix_dim,matrix_dim);

end

% Write the cache record if caching is beneficial
if ismember('op_cache',spin_system.sys.enable)&&(toc>0.1)
    save(filename,'A','-v7.3'); 
end

end

% Input validation function
function grumble(spin_system,operators,spins,operator_type,format)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~(ischar(operators)&&ischar(spins)))&&...
   (~(iscell(operators)&&iscell(spins)))&&...
   (~(ischar(operators)&&isnumeric(spins)))
    error('invalid operator specification.');
end
if iscell(operators)&&iscell(spins)&&(numel(operators)~=numel(spins))
    error('spins and operators cell arrays should have the same number of elements.');
end
if iscell(operators)&&any(~cellfun(@ischar,operators))
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
if ~ischar(operator_type)
    error('operator_type must be a character string.');
end
if ~ischar(format)
    error('format must be a character string.');
end
if ~ismember(format,{'csc','xyz'})
    error('format can be either ''csc'' or ''xyz''.');
end
end

% It is nice to know that the computer understands the 
% problem. But I would like to understand it too.
%
% Eugene Wigner

