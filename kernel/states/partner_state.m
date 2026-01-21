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

% Build partner spin lists and state combinations
nspins=spin_system.comp.nspins;
base_states=repmat({'E'},1,nspins);
base_states{active_spin{2}}=active_spin{1};
partner_spin_lists=cellfun(@(spec) spec{2},partners,'UniformOutput',false);
partner_spins=[partner_spin_lists{:}];
partner_state_sets=cellfun(@(spec) spec{1},partners,'UniformOutput',false);
partner_allowed_state_cells=cellfun(@(states,spins) repmat({states},1,numel(spins)),...
                                    partner_state_sets,partner_spin_lists,...
                                    'UniformOutput',false);
partner_allowed_states=[partner_allowed_state_cells{:}];
% Build partner state combinations
state_counts=cellfun(@numel,partner_allowed_states);
block_sizes=ones(1,numel(state_counts));
if numel(state_counts)>1
    block_sizes(2:end)=cumprod(state_counts(1:end-1));
end
repeat_counts=ones(1,numel(state_counts));
if numel(state_counts)>1
    repeat_counts(1:end-1)=cumprod(state_counts(end:-1:2));
    repeat_counts=repeat_counts(end:-1:1);
end
index_columns=arrayfun(@(count,block,repeat) repmat(kron((1:count).',ones(block,1)),repeat,1),...
                       state_counts,block_sizes,repeat_counts,'UniformOutput',false);
index_matrix=[index_columns{:}];
state_combinations=arrayfun(@(row_idx) build_state_combination(base_states,partner_spins,...
                                        partner_allowed_states,index_matrix(row_idx,:)),...
                            1:size(index_matrix,1),'UniformOutput',false);
spin_list=num2cell(1:nspins);

% Generate the states using Spinach kernel
A=cellfun(@(states) state(spin_system,states,spin_list),state_combinations,'UniformOutput',false);

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
if ~iscell(partners)
    error('partners must be a cell array.');
end
if isempty(partners)
    error('partners must contain at least one specification.');
end
all_partner_spins=[];
for n=1:numel(partners)
    if (~iscell(partners{n}))||(numel(partners{n})~=2)
        error('each partner specification must be a {states,spins} cell array.');
    end
    partner_states=partners{n}{1};
    partner_spin_list=partners{n}{2};
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

% Build the state list for a single partner combination
function current_states=build_state_combination(base_states,partner_spins,...
                                                partner_allowed_states,index_row)
current_states=base_states;
current_states(partner_spins)=cellfun(@(states,idx) states{idx},...
                                      partner_allowed_states,num2cell(index_row),...
                                      'UniformOutput',false);
end

% Challenges improve those who survive.
%
% Frank Herbert
