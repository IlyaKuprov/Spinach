% Overrides default assumptions about interaction terms surviving
% rotating frame transformations. Syntax:
%
%          spin_system=dictum(spin_system,spins,strength)
%
% Parameters:
%
%     spin_system   - Spinach spin system information 
%                     object coming out of assume.m
%
%     spins         - a vector with one or two numbers
%                     or a cell array with one or two
%                     strings, e.g. [2 4] or {'1H'},
%                     where one element would cause
%                     Zeeman interaction assumptions
%                     to be modified, and two elements
%                     would cause coupling assumptions
%                     to be modified.
%
%     strength      - new strength specification, see
%                     the source code of assume.m for
%                     the available strength specs
%
% Outputs:
%
%     spin_system   - updated Spinach spin system in-
%                     formation object that will be
%                     used by hamiltonian.m to build
%                     the Hamiltonian
%
% damien.jeannerat@unige.ch
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dictum.m>

function spin_system=dictum(spin_system,spins,strength)

% Check consistency
grumble(spin_system,spins,strength);

% Override the selected directed coupling records in both directions
if numel(spins)==2
    pairs=spin_system.inter.coupling.strength;
    for q=1:size(pairs,1)
        n=pairs{q,1};
        k=pairs{q,2};
        if isnumeric(spins)
            match=((n==spins(1))&&(k==spins(2)))||((k==spins(1))&&(n==spins(2)));
        else
            isotopes=spin_system.comp.isotopes([n k]);
            match=isequal(isotopes,spins(:)')||isequal(fliplr(isotopes),spins(:)');
        end
        if match
            report(spin_system,sprintf('spins %d, %d: %s changed to %s.',n,k,pairs{q,3},strength));
            pairs{q,3}=strength;
        end
    end
    spin_system.inter.coupling.strength=pairs;

% Numerical specification, Zeeman
elseif isnumeric(spins)&&isscalar(spins)
    
    % Report previous assumption
    report(spin_system,['spin ' num2str(spins) ' (' spin_system.comp.isotopes{spins} ...
                        '), Zeeman assumption: ' spin_system.inter.zeeman.strength{spins}]);
                                 
    % Modify the assumption
    spin_system.inter.zeeman.strength{spins}=strength;
    
    % Report the new assumption
    report(spin_system,['changed on user''s request to: ' spin_system.inter.zeeman.strength{spins}]);
    
% Isotope specification, Zeeman
elseif iscell(spins)&&isscalar(spins)
    
    % Loop over spins pairs
    for n=find(cellfun(@(x)strcmp(spins,x),spin_system.comp.isotopes))
        
        % Report previous assumption
        report(spin_system,['spin ' num2str(n) ' (' spin_system.comp.isotopes{n} ...
                            '), Zeeman assumption: ' spin_system.inter.zeeman.strength{n}]);
        
        % Modify the assumption
        spin_system.inter.zeeman.strength{n}=strength;
        
        % Report the new assumption
        report(spin_system,['changed on user''s request to: ' spin_system.inter.zeeman.strength{n}]);
        
    end
    
else
    
    % Complain and bomb out
    error('incorrect input syntax.');
    
end 

end

% Consistency enforcement
function grumble(spin_system,spins,strength)
if (~isfield(spin_system.inter.coupling,'strength'))||...
   (~isfield(spin_system.inter.zeeman,'strength'))     
    error('assumption information is missing, run assume() before calling this function.');
end
if ~ischar(strength)
    error('strength must be a character string.');
end
if isnumeric(spins)
    if (~isreal(spins))||(~isvector(spins))||(any(spins<1))||...
       (any(mod(spins,1)~=0))||(numel(spins)>2)||(numel(spins)<1)
        error('if numeric, spins must have one or two positive integers.');
    end
    if any(spins>spin_system.comp.nspins)
        error('an element in spins exceeds the number of spins in the system.');
    end
elseif iscell(spins)
    if any(~cellfun(@ischar,spins))
        error('if cell array, spins must contain character strings.');
    end
    if (numel(spins)>2)||(numel(spins)<1)
        error('spins must contain one or two elements.');
    end
    if any(~ismember(spins,spin_system.comp.isotopes))
        error('the requested isotope is not present in the system.');
    end
else
    error('spins must be either a vector of two integers or a cell array of two strings.');
end
end

% It's too bad that stupidity isn't painful.
%
% Anton Szandor LaVey


