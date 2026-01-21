% Partner state expansion; a given state of the specified spin
% is kroneckered with all combinations of the specified states
% of the other specified spins. Syntax:
% 
%       A=partner_state(spin_system,active_spin,partners)
%
% Parameters:
%
%    active_spin - a two-element cell array with the first
%                  element giving the state of the active
%                  spin and the second element giving the
%                  number of the active spin on the main
%                  isotope list, e.g. {'L+',5}
%
%    partners    - a cell array of partner state specifica-
%                  tions of the form
%
%                          { {states_a,spins_a},...
%                            {states_b,spins_b},... }
%
%                  where states_a,b are cell arrays of sta-
%                  tes that the partner spins can have and 
%                  spins_a,b are lists of numbers of those
%                  spins on the main isotope list, e.g.
%
%                         { {{'E'  'Lz'},[1 5]},...
%                           {{'L+' 'L-'},[2 7]},... }
%
% Outputs:
%
%    A - a cell array of spin states (matrices in Hilbert space,
%        vectors in Liouville space) with the active spin in the
%        specified state and the partner spins in all combinati-
%        ons specified by the user. All spins not explicitly
%        mentioned in the input will be in their 'E' states.
%
% Example: in a five-spin system, the following call
%
%  A=partner_state(spin_system,{'L+',2},{{{'E','Lz'},[1 3]}})
%
%          will return the following state array
%
%  A={state(spin_system,{'E' ,'L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'Lz','L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'E' ,'L+','Lz','E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'Lz','L+','Lz','E' ,'E'},{1 2 3 4 5})};
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=partner_state.m>

function A=partner_state(spin_system,active_spin,partners)

% Check consistency
grumble(spin_system,active_spin,partners);

% Normalize partner specifications
if isempty(partners)
    partner_specs={};
elseif (iscell(partners))&&(numel(partners)==2)&&...
       (iscell(partners{1}))&&(isnumeric(partners{2}))
    partner_specs={partners};
else
    partner_specs=partners;
end

% Build partner spin lists and state combinations
nspins=spin_system.comp.nspins;
base_states=repmat({'E'},1,nspins);
base_states{active_spin{2}}=active_spin{1};
partner_spins=[];
partner_allowed_states={};
for n=1:numel(partner_specs)
    partner_states=partner_specs{n}{1};
    partner_spin_list=partner_specs{n}{2};
    for k=1:numel(partner_spin_list)
        partner_spins=[partner_spins partner_spin_list(k)];
        partner_allowed_states{end+1}=partner_states;
    end
end
if isempty(partner_spins)
    index_matrix=zeros(1,0);
else
    state_counts=cellfun(@numel,partner_allowed_states);
    index_axes=cell(1,numel(state_counts));
    for n=1:numel(state_counts)
        index_axes{n}=1:state_counts(n);
    end
    [index_axes{:}]=ndgrid(index_axes{:});
    index_matrix=zeros(numel(index_axes{1}),numel(state_counts));
    for n=1:numel(state_counts)
        index_matrix(:,n)=index_axes{n}(:);
    end
end
state_combinations=cell(size(index_matrix,1),1);
for n=1:size(index_matrix,1)
    current_states=base_states;
    for k=1:numel(partner_spins)
        current_states{partner_spins(k)}=partner_allowed_states{k}{index_matrix(n,k)};
    end
    state_combinations{n}=current_states;
end
spin_list=num2cell(1:nspins);

% Generate the states using Spinach kernel
A=cell(numel(state_combinations),1);
for n=1:numel(state_combinations)
    A{n}=state(spin_system,state_combinations{n},spin_list);
end

end

% Consistency enforcement
function grumble(spin_system,active_spin,partners)
if (~iscell(active_spin))||(numel(active_spin)~=2)||...
   (~ischar(active_spin{1}))||(~isnumeric(active_spin{2}))||...
   (~isscalar(active_spin{2}))||(~isreal(active_spin{2}))||...
   (mod(active_spin{2},1)~=0)||(active_spin{2}<1)
    error('active_spin must be a {state,spin_number} cell array.');
end
if active_spin{2}>spin_system.comp.nspins
    error('active spin index exceeds the number of spins.');
end
if isempty(partners)
    return
end
if ~iscell(partners)
    error('partners must be a cell array.');
end
if (numel(partners)==2)&&(iscell(partners{1}))&&(isnumeric(partners{2}))
    partner_specs={partners};
else
    partner_specs=partners;
end
all_partner_spins=[];
for n=1:numel(partner_specs)
    if (~iscell(partner_specs{n}))||(numel(partner_specs{n})~=2)
        error('each partner specification must be a {states,spins} cell array.');
    end
    partner_states=partner_specs{n}{1};
    partner_spin_list=partner_specs{n}{2};
    if (~iscell(partner_states))||(~all(cellfun(@ischar,partner_states)))
        error('partner states must be a cell array of character strings.');
    end
    if (~isnumeric(partner_spin_list))||(~isvector(partner_spin_list))||...
       any(mod(partner_spin_list,1)~=0,'all')||any(partner_spin_list<1,'all')
        error('partner spin lists must contain positive integers.');
    end
    all_partner_spins=[all_partner_spins partner_spin_list(:).'];
end
if numel(unique(all_partner_spins))~=numel(all_partner_spins)
    error('partner spin lists must not contain duplicates.');
end
if any(all_partner_spins>spin_system.comp.nspins,'all')
    error('partner spin index exceeds the number of spins.');
end
if any(all_partner_spins==active_spin{2},'all')
    error('active spin must not appear in partner spin lists.');
end
end

% Challenges improve those who survive.
%
% Frank Herbert
