% Converts user-friendly descriptions of spin states and operators into the 
% formal description (opspec) used by Spinach kernel. The function supports
% three types of syntax:
%
% 1. If both inputs are strings
%
%       [opspecs,coeffs]=human2opspec(spin_system,'Lz','13C')
%
% the function returns a list of single-spin opspecs for all spins with the
% specified name. In the example above, the list of Lz operator specificati-
% ons for all 13C nuclei in the system would be returned. Valid labels for 
% operators in this type of call are 'E' (identity), 'Lz', 'Lx', 'Ly', 'L+',
% 'L-', and 'Tl,m' (irreducible spherical tensor, l and m are integers). Va-
% lid labels for the spins are standard isotope names as well as 'electrons', 
% 'nuclei', and 'all'.
%
% 2. If one input is a string and the other is a vector
%
%       [opspecs,coeffs]=human2opspec(spin_system,'Lz',[1 2 4])
%
% the function returns a list of single-spin opspecs for all spins with the
% specified number. In the example above, the list of Lz operator specifica-
% tions for all 13C nuclei in the system would be returned. Valid labels for
% operators are the same as in Item 1 above.
%
% 3. If the two inputs are a cell array of strings and a cell array of num-
%    bers, a product operator specification is produced
% 
%       [opspecs,coeffs]=human2opspec(spin_system,{'Lz','Ly'},{1,2})
%
% would return the Lz(x)Ly product operator specification with Lz on spin 1
% and Ly on spin 2. Valid labels for operators are the same as in Item 1.
%
% Outputs:
%
%    opspecs    - Spinach operator specification: a cell array of
%                 row vectors specifying which operator enters the
%                 Kronecker product for which spin.
%
%    coeffs     - coefficient with which each of the Kronecker pro-
%                 ducts enters the overall sum.
%
% Notes: direct calls to this function are not necessary, use operator.m and
%        state.m functions instead.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=human2opspec.m>

function [opspecs,coeffs]=human2opspec(spin_system,operators,spins)

% Check input validity
grumble(operators,spins);

% Preliminary operator inspection
if ischar(operators)&&strcmp(operators,'Lx')

    % Recursive call for Lx = (Lp + Lm)/2
    [opspecs_a,coeffs_a]=human2opspec(spin_system,'L+',spins);
    [opspecs_b,coeffs_b]=human2opspec(spin_system,'L-',spins);
    opspecs=[opspecs_a; opspecs_b]; coeffs=[coeffs_a; coeffs_b]/2; return

elseif ischar(operators)&&strcmp(operators,'Ly')

    % Recursive call for Ly = (Lp - Lm)/2i
    [opspecs_a,coeffs_a]=human2opspec(spin_system,'L+',spins);
    [opspecs_b,coeffs_b]=human2opspec(spin_system,'L-',spins);
    opspecs=[opspecs_a; opspecs_b]; coeffs=[coeffs_a; -coeffs_b]/2i; return

% Likewise for product states
elseif iscell(operators)

    % Inspect each term
    for n=1:numel(operators)

        if strcmp(operators{n},'Lx')

            % Recursive call for Lx = (Lp + Lm)/2
            operators{n}='L+'; [opspecs_a,coeffs_a]=human2opspec(spin_system,operators,spins);
            operators{n}='L-'; [opspecs_b,coeffs_b]=human2opspec(spin_system,operators,spins);
            opspecs=[opspecs_a; opspecs_b]; coeffs=[coeffs_a; coeffs_b]/2; return

        elseif strcmp(operators{n},'Ly')

            % Recursive call for Ly = (Lp - Lm)/2i
            operators{n}='L+'; [opspecs_a,coeffs_a]=human2opspec(spin_system,operators,spins);
            operators{n}='L-'; [opspecs_b,coeffs_b]=human2opspec(spin_system,operators,spins);
            opspecs=[opspecs_a; opspecs_b]; coeffs=[coeffs_a; -coeffs_b]/2i; return

        end

    end

end

% 'Lz','1H' type call returns a sum
if ischar(operators)&&ischar(spins)
    
    % Parse the specification
    switch spins
        
        case 'all'
            
            % Use all spins in the system
            spin_numbers=1:spin_system.comp.nspins;
            
        case 'electrons'
            
            % Include electrons of any multiplicity
            spin_numbers=find(cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
            
        case 'nuclei'
            
            % Include non-electrons
            spin_numbers=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
            
        otherwise
            
            % Count specific spins
            spin_numbers=find(strcmp(spin_system.comp.isotopes,spins));
            
    end
    
    % Bomb out if the list of spins is empty
    if numel(spin_numbers)==0, error('no such spins in the system.'); end
    
    % Preallocate descriptors
    coeffs=zeros(numel(spin_numbers),1);
    opspecs=cell(numel(spin_numbers),1);
    
    % Use recursive calls
    for n=1:numel(spin_numbers)
        [opspecs(n),coeffs(n)]=human2opspec(spin_system,{operators},{spin_numbers(n)});
    end

    % Return the outcome
    return

% 'Lz',[1 2 4] type call returns a sum
elseif ischar(operators)&&isnumeric(spins)
    
    % Preallocate descriptors
    coeffs=zeros(numel(spins),1);
    opspecs=cell(numel(spins),1);
    
    % Use recursive calls
    for n=1:numel(spins)
        [opspecs(n),coeffs(n)]=human2opspec(spin_system,{operators},{spins(n)});
    end
    
    % Return the outcome
    return
    
% {'Lz','L+'},{1,4} type call returns an operator
elseif iscell(operators)&&iscell(spins)

    % Start with an empty opspec and a unit coefficient
    opspecs={zeros(1,spin_system.comp.nspins)}; coeffs=1;
    
    % Parse operator selection
    for n=1:numel(operators)
        
        % Operator type
        switch operators{n}

            case 'Em'

                % Empty cavity state
                opspecs{1}(spins{n})=-4;

            case 'An'
                
                % Annihilation operator
                opspecs{1}(spins{n})=-3;

            case 'Nu'
                
                % Number operator
                opspecs{1}(spins{n})=-2;

            case 'Cr'
                
                % Creation operator
                opspecs{1}(spins{n})=-1;

            case 'E'
                
                % Unit operator
                opspecs{1}(spins{n})=0;
                
            case 'L+'
                
                % Raising operator
                opspecs{1}(spins{n})=1;

                % T(1,+1) coefficient
                coeffs=-sqrt(2)*coeffs;
                
            case 'Lz'
                
                % Z projection operator
                opspecs{1}(spins{n})=2;
                
            case 'L-'
                
                % Lowering operator
                opspecs{1}(spins{n})=3;
                
                % T(1,-1) coefficient
                coeffs=sqrt(2)*coeffs;

            otherwise
                
                % Validate the irreducible spherical tensor input
                if isempty(regexp(operators{n},'^T([\+\-]?\d+),([\+\-]?\d+)$','once'))
                    error('unrecognized operator specification.');
                end
                
                % Extract the quantum numbers
                indices=textscan(operators{n},'T%n,%n'); l=indices{1}; m=indices{2};
                
                % Validate the quantum numbers
                if (l<0)||(abs(m)>l)
                    error('invalid indices in irreducible spherical tensors.');
                end
                
                % Write the specification
                opspecs{1}(spins{n})=lm2lin(l,m); coeffs=1*coeffs;
                
        end
        
    end
    
else
    
    % Complain and bomb out
    error('invalid operator or state specification.');
    
end

end

% Input validation
function grumble(operators,spins)
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
end

% A "moral commandment" is a contradiction in terms. The moral is the
% chosen, not the forced; the understood, not the obeyed. The moral is
% the rational, and reason accepts no commandments.
%
% Ayn Rand

