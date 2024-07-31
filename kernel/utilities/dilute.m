% Splits the spin system into several independent subsystems, each
% containing only one instance of a user specified isotope or iso-
% tope tuple that is deemed "dilute". All spin system data is upda-
% ted accordingly. Basis set information, if found, is destroyed.
% Syntax:
%
%          subsystems=dilute(spin_system,isotope,tuples)
%
% Parameters:
%
%    spin_system   - spin system object obtained 
%                    from create.m function
%
%    isotope       - isotope specification string,
%                    for example '13C'
%
%    tuples        - 1 returns all subsystems with
%                    a single instance of the isoto-
%                    pe, 2 returns all subsystems
%                    where two instances of the iso-
%                    tope are present, etc.
%
% Outputs:
%
%    subsystems    - a cell array of spin system 
%                    objects, with each cell cor-
%                    responding to one of the new-
%                    ly formed isotopomers.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
% bud.macaulay@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dilute.m>

function subsystems=dilute(spin_system,isotope,tuples)

% Default is singles
if ~exist('tuples','var'), tuples=1; end

% Check consistency
grumble(isotope,tuples);

% Inform the user
report(spin_system,['treating ' isotope ' as a dilute, ' ... 
                    'picking tuples of size ' num2str(tuples)]);

% Find out which spins belong to the dilute species
dilute_spins=find(cellfun(@(x)strcmp(x,isotope),spin_system.comp.isotopes));

% Report the quantity
report(spin_system,[num2str(numel(dilute_spins)) ' instances of '...
                    isotope ' found in the system.']);
       
% Disallow picking more than there is
if numel(dilute_spins)<tuples
    error(['not enough ' isotope ' in the system to make the selection.']);
end

% Get the combinations
combos=nchoosek(dilute_spins,tuples);
combos=unique(sort(combos,2),'rows');

% Preallocate the answer
subsystems=cell(size(combos,1),1);

% Create new spin systems
for n=1:size(combos,1)
    subsystems{n}=kill_spin(spin_system,setdiff(dilute_spins,combos(n,:)));
end

% Inform the user
report(spin_system,[num2str(numel(subsystems)) ' spin subsystems returned.']); 

end

% Consistency enforcement
function grumble(isotope,tuples)
if ~ischar(isotope)
    error('isotope must be a character string.');
end
if (~isnumeric(tuples))||(~isreal(tuples))||...
   (~isscalar(tuples))||(tuples<1)||(mod(tuples,1)~=0)
    error('tuples must be a positive integer.');
end
end

% "Dear Mother, good news today."
%
% Albert Einstein, in a 1919 postcard to his mother telling her that his
% general theory of relativity had been proven.

