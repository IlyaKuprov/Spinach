% Test of internal consistency for state and operator 
% generation across the three formalisms supported by
% Spinach. Two- and four-spin states are tested.
%
% i.kuprov@soton.ac.uk

function state_consistency()

% Magnet field
sys.magnet=14.1;

% Set the spin system
sys.isotopes={'1H','1H','1H','1H'};
inter.zeeman.scalar={0 0 0 0};

% Loop over the formalisms
Fs={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
for n=1:numel(Fs)

    % Basis set
    bas.formalism=Fs{n};
    bas.approximation='none';

    % Hush the logs
    sys.output='hush';
    sys.disable={'hygiene'};

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Unit state
    U=state(spin_system,{'E','E','E','E'},{1 2 3 4});

    % Two-spin singlet-triplet state sum test
    combos={[1 2],[1 3],[1 4],[2 1],[2 3],[2 4],...
            [3 1],[3 2],[3 4],[4 1],[4 2],[4 3]};
    parfor k=1:numel(combos) %#ok<*PFBNS>
        S=singlet(spin_system,combos{k}(1),combos{k}(2));
        [TU,T0,TD]=triplet(spin_system,combos{k}(1),combos{k}(2));
        if norm(S+TU+T0+TD-U,1)>1e-6
            error([Fs{n} ': state construction test FAILED.']);
        end
    end
    
    % Four-spin singlet-triplet state sum test
    combos=num2cell(perms([1 2 3 4]),2);
    parfor k=1:numel(combos)
        S=four_spin_states(spin_system,num2cell(combos{k}),'S(x)S')+...
          four_spin_states(spin_system,num2cell(combos{k}),'S(x)TU')+...
          four_spin_states(spin_system,num2cell(combos{k}),'S(x)T0')+...
          four_spin_states(spin_system,num2cell(combos{k}),'S(x)TD')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TU(x)S')+...
          four_spin_states(spin_system,num2cell(combos{k}),'T0(x)S')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TD(x)S')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TU(x)TU')+...
          four_spin_states(spin_system,num2cell(combos{k}),'T0(x)TU')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TD(x)TU')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TU(x)T0')+...
          four_spin_states(spin_system,num2cell(combos{k}),'T0(x)T0')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TD(x)T0')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TU(x)TD')+...
          four_spin_states(spin_system,num2cell(combos{k}),'T0(x)TD')+...
          four_spin_states(spin_system,num2cell(combos{k}),'TD(x)TD');
        if norm(S-U,1)>1e-6
            error([Fs{n} ': state construction test FAILED.']); 
        end
    end

    % Good news to the user
    disp([Fs{n} ': state construction test PASSED.']);

end

end

