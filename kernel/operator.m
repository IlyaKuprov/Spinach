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
% Note: operator caching is supported, add 'op_cache' to sys.enable array 
%       to enable; make sure your scratch storage is fast.
%
% ledwards@cbs.mpg.de
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=operator.m>

function A=operator(spin_system,operators,spins,operator_type,format)

% The default superoperator type is commutation
if ~exist('operator_type','var'), operator_type='comm'; end

% The default sparse matrix format is CSC
if ~exist('format','var'), format='csc'; end

% Check consistency and start the clock
grumble(spin_system,operators,spins,...
        operator_type,format); tic;

% Load the cache record if one exists
if ismember('op_cache',spin_system.sys.enable)

    % Combine specification, isotopes, and basis hash
    op_hash=md5_hash({operators,spins,operator_type,...
                      format,spin_system.comp.iso_hash,...
                      spin_system.bas.basis_hash});

    % Generate the cache record name in the scratch directory
    filename=[spin_system.sys.scratch filesep 'spinach_op_' op_hash '.mat'];

    % Check if the file exists
    if exist(filename,'file')

        % Try to use
        try

            % Try to load
            load(filename,'A');

            % Check load success
            if exist('A','var')
                return; 
            else
                % Do not make a fuss on fail
            end

        catch

            % Do not make a fuss on fail
            
        end

    end
    
end

% Parse the human specification into Spinach notation
[opspecs,coeffs]=human2opspec(spin_system,operators,spins);

% Start a cell array
A=cell(numel(opspecs),1);

% This is formalism-dependent
switch spin_system.bas.formalism
    
    % Spherical tensors
    case 'sphten-liouv'

        % Parallelisation efficiency
        types=spin_system.comp.types;
        
        % Build summation terms
        parfor n=1:numel(opspecs)

            % Check physics type
            if all(strcmp('S',types),'all')

                % For spins: get the superoperator in XYZ format
                A{n}=p_superop(spin_system,opspecs{n},operator_type);

                % For spins: apply coefficient
                A{n}(:,3)=coeffs(n)*A{n}(:,3);

            else

                % For cavities and phonons: Sarbojoy's PhD project :) :) :)
                error('Cavities and phonons not yet available in sphten-liouv.');

            end

        end
                
    % Zeeman basis formalisms
    case {'zeeman-wavef',...
          'zeeman-hilb',...
          'zeeman-liouv'}

       
        % Parallelisation efficiency
        mults=spin_system.comp.mults;
        formalism=spin_system.bas.formalism;
        types=spin_system.comp.types;

        % Build summation terms
        parfor n=1:numel(opspecs)

            % Index pertinent spins
            active_spins=find(opspecs{n});

            % Find out which substance hosts the opspec
            % and narrow it down to that substance
            subst=which_subst(spin_system,active_spins);
            spins_in_subst=spin_system.chem.parts{subst}(:);
            opspec_in_subst=opspecs{n}(spins_in_subst);

            % Pull out spin multiplicities in the substance
            mults_in_subst=spin_system.comp.mults(spins_in_subst);

            % Start the kron
            B=sparse(coeffs(n));
            
            % Over the substance content
            for k=1:numel(opspec_in_subst)

                % Check physics type
                switch types{k} %#ok<PFBNS>
                
                    case 'S' % Spins

                    % Spin: get spin state indices
                    [L,M]=lin2lm(opspec_in_subst(k));

                    % Spin: get irreducible spherical tensors
                    ist=irr_sph_ten(mults_in_subst(k),L);

                        % Update the kron
                        B=kron(B,ist{L-M+1});

                    case {'C','V','T'} % Cavities, phonons, transmons

                    % Cavities and phonons: get operator
                    T=boson_oper(mults_in_subst(k),...
                                 opspec_in_subst(k));

                        % Update the kron
                        B=kron(B,bmon{opspecs{n}(k)+1});

                    otherwise

                        % Complain and bomb out
                        error('unknown particle type.');

                end
                
            end
            
            % Move to Liouville space if necessary
            if strcmp(formalism,'zeeman-liouv')
                B=hilb2liouv(B,operator_type);
            end

            % Convert sparse array from CSC to XYZ indexing
            [rows,cols,vals]=find(B); A{n}=[rows cols vals];

            % Get the global multi-substance basis index offset
            idx_offset=sum(spin_system.bas.nstates(1:(subst-1)));

            % Apply the offset
            A{n}(:,1)=A{n}(:,1)+idx_offset;
            A{n}(:,2)=A{n}(:,2)+idx_offset;

        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism.');
        
end

% XYZ format sum
A=cell2mat(A);

% Convert to CSC format
if strcmp(format,'csc')

    % Get the global operator dimension
    matrix_dim=sum(spin_system.bas.nstates);

    % Make a sparse matrix, making sure it is complex for later
    A=sparse(A(:,1),A(:,2),complex(A(:,3)),matrix_dim,matrix_dim);

end

% Write the cache record if caching is beneficial
if ismember('op_cache',spin_system.sys.enable)&&(toc>0.1)

    % Do not fight other workers
    if ~exist(filename,'file')

        % Try to save
        try

            % Modern format, compressed
            save(filename,'A','-v7.3'); drawnow;

        catch

            % Do not make a fuss on fail, this can happen
            % for large parallel pools where many workers
            % may be trying to write the same file.

        end

    end

end

end

% Consistency enforcement
function grumble(spin_system,operators,spins,operator_type,format)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() first.');
end
if (~(ischar(operators)&&ischar(spins)))&&...
   (~(iscell(operators)&&iscell(spins)))&&...
   (~(ischar(operators)&&isnumeric(spins)))
    error('invalid operator specification.');
end
if iscell(operators)&&iscell(spins)&&(numel(operators)~=numel(spins))
    error('spins and operators cell arrays must have the same number of elements.');
end
if iscell(operators)&&any(~cellfun(@ischar,operators))
    error('all elements of the operators cell array must be strings.');
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

