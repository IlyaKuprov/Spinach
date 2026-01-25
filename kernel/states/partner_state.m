% Partner state expansion; a given state of the specified spins
% is kroneckered with all combinations of the specified states
% of specfied partner spins. Syntax:
% 
%  [A,descr]=partner_state(spin_system,active_spin,partners)
%
% Parameters:
%
%    set_spin    - a cell array of two-element cell arrays 
%                  with the first element giving the state
%                  of the spin and the second element num-
%                  ber on the isotope list, for example
%
%                                {{'L+',3}}
%
%                  These spins will have their states set
%                  immutably as specified.
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
%                  These spins will have their state varied
%                  combinatorially as specified.
%
% Outputs:
%
%    A - a cell array of spin states (matrices in Hilbert space,
%        vectors in Liouville space) with the active spin in the
%        specified state and the partner spins in all combinati-
%        ons specified by the user. All spins not explicitly
%        mentioned in the input will be in their 'E' states.
%
%    descr - a cell array of product structure descriptors for 
%            each element of A
%
% Example: in a five-spin system, the following call
%
%  A=partner_state(spin_system,{{'L+',2}},{{{'E','Lz'},[1 3]}})
%
%          will return the following state array
%
%  A={state(spin_system,{'E' ,'L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'Lz','L+','E' ,'E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'E' ,'L+','Lz','E' ,'E'},{1 2 3 4 5}),...
%     state(spin_system,{'Lz','L+','Lz','E' ,'E'},{1 2 3 4 5})};
%
%          and the following descriptor array
%
%     descr={{'E' ,'L+','E' ,'E' ,'E'},...
%            {'Lz','L+','E' ,'E' ,'E'},...
%            {'E' ,'L+','Lz','E' ,'E'},...
%            {'Lz','L+','Lz','E' ,'E'}};
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=partner_state.m>

function [A,descr]=partner_state(spin_system,set_spin,partners)

% Write the states of the spins that stay unchanged
base_state=repmat({'E'},1,spin_system.comp.nspins);
for n=1:numel(set_spin)
    base_state{set_spin{n}{2}}=set_spin{n}{1};
end

% Pull out partner spin lists and partner state sets
partner_spin_lists=cellfun(@(spec)spec{2},partners,'UniformOutput',false);
partner_state_sets=cellfun(@(spec)spec{1},partners,'UniformOutput',false);

% Compile the active partner list
active_partners=[partner_spin_lists{:}];

% Compile state lists for each partner
partner_state_sets=cellfun(@(x,y)repmat({x},1,numel(y)),...
                           partner_state_sets,...
                           partner_spin_lists,'UniformOutput',false);
partner_state_sets=[partner_state_sets{:}];

% Start descriptor
descr={base_state};

% Over active partner spins
for n=1:numel(active_partners)

    % Pull current partner
    partner=active_partners(n);
    
    % Over the partner state set
    for k=1:numel(partner_state_sets{n})

        % Pull current partner spin state
        partner_state=partner_state_sets{n}{k};

        % Scan descriptors
        for m=1:numel(descr)

            % Pull descriptor line
            current_line=descr{m};

            % Modify the descriptor line if necessary
            if ~strcmp(current_line{partner},partner_state)

                % Write in the new state and append
                current_line{partner}=partner_state;
                descr=[descr; {current_line}]; %#ok<AGROW>

            end

        end

    end

end

% Get the full spin list    
full_spin_list=num2cell(1:spin_system.comp.nspins);

% Generate the states
A=cell(1,numel(descr));
parfor n=1:numel(descr)
    A{n}=state(spin_system,descr{n},full_spin_list);
end

end

% Challenges improve those who survive.
%
% Frank Herbert

