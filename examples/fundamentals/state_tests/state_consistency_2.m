% Test of consistency in the projection between spherical
% tensor basis set and Zeeman basis set.
%
% ilya.kuprov@weizmann.ac.il

function state_consistency_2()

% Magneti field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'14N','235U'};

% No interactions
inter.zeeman.scalar={0 0};

% Hilbert space, Zeeman basis
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=create(sys, inter);
spin_system=basis(spin_system, bas);

% A suitably complicated state
Op{1}=state(spin_system,{'Lz','Lx'},{1,2})+...
      state(spin_system,{'L+'},{1});

% Liouville space, Zeeman basis
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=create(sys, inter);
spin_system=basis(spin_system, bas);

% A suitably complicated state
Op{2}=state(spin_system,{'Lz','Lx'},{1,2})+...
      state(spin_system,{'L+'},{1});

% Fold back into Hilbert space
Op{2}=reshape(Op{2},[24 24]);

% Liouville space, IST basis
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys, inter);
spin_system=basis(spin_system, bas);

% A suitably complicated state
Op{3}=state(spin_system,{'Lz','Lx'},{1,2})+...
      state(spin_system,{'L+'},{1});

% Project into Zeeman basis
Op{3}=sphten2zeeman(spin_system)*Op{3};

% Fold back into Hilbert space
Op{3}=reshape(Op{3},[24 24]);

% Check the results
if (norm(Op{1}-Op{2},1)>1e-6)||...
   (norm(Op{2}-Op{3},1)>1e-6)||...
   (norm(Op{3}-Op{1},1)>1e-6)

    % Complain and bomb out
    error('State consistency test FAILED.');

else

    % Good news to the user
    disp('State consistency test PASSED.');

end

end

