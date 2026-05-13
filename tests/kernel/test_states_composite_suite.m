% Tests composite state generators in kernel/states. Syntax:
%
%                    result=test_states_composite_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks unit-state normalisation, two-spin singlet/triplet
% projectors, four-spin product states, and partner-state enumeration.
%
% ilya.kuprov@weizmann.ac.il

function result=test_states_composite_suite()

% Announce the test target
fprintf('TESTING: Composite state generators\n');

% State the state-generation target of the test
result=new_test_result('kernel/states_composite_suite',...
                       'Composite state generators',...
                       'state helper functions must produce the textbook density operators.');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Check the Hilbert-space thermodynamic unit state
result=test_close(result,'unit_state zeeman-hilb',unit_state(spin_system),speye(2),1e-15,1e-15,...
                  'the Hilbert-space unit state is the identity density matrix');

% Check the Zeeman-Liouville unit vector normalisation
bas.formalism='zeeman-liouv';
spin_system=test_spin_system(sys,inter,bas);
unit=speye(2); unit=unit(:); unit=unit/norm(unit,2);
result=test_close(result,'unit_state zeeman-liouv',unit_state(spin_system),unit,1e-15,1e-15,...
                  'the Zeeman-Liouville unit state is the normalised vectorised identity');

% Build a two-proton Hilbert-space spin system
clear sys inter bas
sys.magnet=0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
inter.coupling.scalar{1,2}=0;
inter.coupling.scalar{2,2}=0;
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Textbook two-spin wavefunctions in the Zeeman product basis
alpha=[1;0];
beta=[0;1];
aa=kron(alpha,alpha);
ab=kron(alpha,beta);
ba=kron(beta,alpha);
bb=kron(beta,beta);
sing=(ab-ba)/sqrt(2);
trip_zero=(ab+ba)/sqrt(2);

% Check singlet and triplet projectors
result=test_close(result,'singlet projector',singlet(spin_system,1,2),sing*sing',1e-14,1e-14,...
                  'the singlet state is |alpha beta-beta alpha><...|/2');
[TU,T0,TD]=triplet(spin_system,1,2);
result=test_close(result,'triplet-up projector',TU,aa*aa',1e-14,1e-14,...
                  'the triplet-up state is |alpha alpha><alpha alpha|');
result=test_close(result,'triplet-zero projector',T0,trip_zero*trip_zero',1e-14,1e-14,...
                  'the triplet-zero state is |alpha beta+beta alpha><...|/2');
result=test_close(result,'triplet-down projector',TD,bb*bb',1e-14,1e-14,...
                  'the triplet-down state is |beta beta><beta beta|');

% Build a four-proton Hilbert-space spin system
clear sys inter bas
sys.magnet=0;
sys.isotopes={'1H','1H','1H','1H'};
inter.zeeman.scalar={0,0,0,0};
inter.coupling.scalar{1,2}=0;
inter.coupling.scalar{2,3}=0;
inter.coupling.scalar{3,4}=0;
inter.coupling.scalar{4,4}=0;
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Check the product of two singlet projectors
singlet_pair=sing*sing';
result=test_close(result,'four_spin_states S(x)S',four_spin_states(spin_system,[1 2 3 4],'S(x)S'),...
                  kron(singlet_pair,singlet_pair),1e-13,1e-13,...
                  'S(x)S is the tensor product of singlets on spins 1-2 and 3-4');

% Build a three-proton Hilbert-space spin system for partner-state enumeration
clear sys inter bas
sys.magnet=0;
sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={0,0,0};
inter.coupling.scalar{1,2}=0;
inter.coupling.scalar{2,3}=0;
inter.coupling.scalar{3,3}=0;
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Enumerate two partner spins that may each be E or Lz while spin 2 is L+
[A,descr]=partner_state(spin_system,{{'L+',2}},{{{'E','Lz'},[1 3]}});
expected_descr={{'E','L+','E'};...
                {'Lz','L+','E'};...
                {'E','L+','Lz'};...
                {'Lz','L+','Lz'}};
result=test_true(result,'partner_state descriptor order',isequal(descr(:),expected_descr(:)),...
                 'partner_state() enumerates the Cartesian product of allowed partner states');

% Check every enumerated partner state against a direct state() call
full_spin_list={1,2,3};
for n=1:numel(expected_descr)
    rho_ref=state(spin_system,expected_descr{n},full_spin_list);
    result=test_close(result,['partner_state element ' int2str(n)],A{n},rho_ref,1e-14,1e-14,...
                      'each partner-state descriptor must map to the corresponding direct product state');
end

end


