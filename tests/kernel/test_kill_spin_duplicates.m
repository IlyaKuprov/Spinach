% Tests rejection of repeated particle indices. Syntax:
%
%                    result=test_kill_spin_duplicates()
%
% Outputs:
%
%     result - duplicate-index validation test results
%
% Removing particle two from a two-proton system must leave the same
% one-proton system as an independent create() call. Repeating that index
% must fail before changing the particle count or chemical partitioning.
% Exact records and roundoff-level Hamiltonian agreement are required.
%
% ilya.kuprov@weizmann.ac.il

function result=test_kill_spin_duplicates()

% Announce the particle-removal contract
fprintf('TESTING: repeated particle indices\n');
result=new_test_result('kernel/kill_spin_duplicates','Repeated particle indices',...
                      'Repeated indices cannot corrupt the particle count.');

% Construct the original and independent retained systems
sys.magnet=1; sys.isotopes={'1H','1H'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter=struct(); bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=create(sys,inter); rejected=false;
try
    kill_spin(spin_system,[2 2]);
catch err
    rejected=contains(err.message,'hit_list')&&contains(err.message,'repeated');
end
result=test_true(result,'duplicate indices rejected',rejected,...
                 'The validation error must identify repeated hit_list indices.');

% Unique indices and logical masks retain their original semantics
sys.isotopes={'1H'}; reference=assume(basis(create(sys,inter),bas),'labframe');
for hit_list={2,[false true]}
    observed=assume(basis(kill_spin(spin_system,hit_list{1}),bas),'labframe');
    result=test_true(result,'retained particle records',...
                     isequal(observed.comp.isotopes,reference.comp.isotopes)&&...
                     observed.comp.nspins==reference.comp.nspins,...
                     'A unique removal leaves exactly the independent one-proton system.');
    result=test_true(result,'retained chemical records',...
                     isequaln(observed.chem,reference.chem),...
                     'The retained chemical partition contains particle one exactly once.');
    result=test_close(result,'retained Hamiltonian',...
                      hamiltonian(observed),hamiltonian(reference),...
                      100*eps(),100*eps(),'Removing a spectator preserves the one-spin Hamiltonian.');
end

end


