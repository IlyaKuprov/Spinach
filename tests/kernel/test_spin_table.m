% Tests sourced isotope metadata and numerical simulation boundaries. Syntax:
%
%                         result=test_spin_table()
%
% Outputs:
%
%    result - source conventions, state identity, particle signs, and API checks
%
% talos@spindynamics.org

function result=test_spin_table()

% Obtain the physical table through the public API
result=new_test_result('kernel/spin_table','Literature isotope table',...
                       'metadata must retain source qualifications without guessing physical inputs.');
[~,~,data]=spin('table');
result=test_true(result,'unique physical keys',...
                 numel(unique(data.isotope))==height(data),...
                 'each observed nuclear state or particle has exactly one key');
result=test_true(result,'academic row attribution',all(contains(data.source_dois,"10.")),...
                 'every physical row names an academic DOI source');
result=test_true(result,'required recent nuclei',...
                 all(ismember({'264Lr','268Sg','231Am'},data.Properties.RowNames)),...
                 'published later discoveries with measured lifetimes above one second are retained');

% Compare the shipped payload with an independently rebuilt temporary source copy
temp_dir=tempname(); mkdir(temp_dir);
cleanup=onCleanup(@()rmdir(temp_dir,'s'));
source_file=fullfile(fileparts(which('spin')),'..','etc','isotopes.tsv');
copyfile(source_file,fullfile(temp_dir,'isotopes.tsv'));
rebuilt=build_isotopes(fullfile(temp_dir,'isotopes.tsv'));
result=test_true(result,'release source equality',isequaln(data,rebuilt),...
                 'the shipped MAT must represent every TSV value and source qualification');
clear cleanup

% Check exact spin-zero restrictions and distinct tantalum nuclear states
[gamma,mult,zero]=spin('12C');
result=test_true(result,'spin zero',gamma==0&&mult==1&&zero.quad_moment==0,...
                 'spin zero has no magnetic or quadrupole moment');
ground=data('180Ta',:); isomer=data('180Ta_m',:);
result=test_true(result,'tantalum state identity',...
                 ground.spin==1&&isomer.spin==9&&ground.half_life<86400&&...
                 isinf(isomer.half_life)&&isomer.abundance>0&&isnan(ground.abundance),...
                 'natural tantalum-180 abundance belongs to the long-lived isomer');

% Preserve interval abundance rather than inventing a midpoint
hydrogen=data('1H',:);
result=test_true(result,'hydrogen interval',...
                 isnan(hydrogen.abundance)&&hydrogen.abund_low==0.99972&&...
                 hydrogen.abund_high==0.99999,...
                 'the adopted terrestrial hydrogen composition is a published interval');

% Check particle charge signs and the mean-lifetime conversion
[elec,~]=spin('E'); [posi,~]=spin('E+');
[muon,~,negative]=spin('M'); [anti_muon,~]=spin('M+');
result=test_true(result,'charged counterparts',elec<0&&posi==-elec&&...
                 muon<0&&anti_muon==-muon,...
                 'electron and negative muon signs differ from their CPT counterparts');
result=test_close(result,'free muon half-life',negative.half_life,...
                  log(2)*2.1969811e-6,1e-20,1e-12,...
                  'the source mean lifetime is not labelled as a half-life');
result=test_true(result,'particle classification',...
                 ~isnucleus('E+')&&~isnucleus('M+')&&~isnucleus('anti1H')&&...
                 isnucleus('13C'),...
                 'physical particle rows must not be mistaken for nuclei');

% Anchor the spin-three-halves hyperon conversion and particle g-factor convention
[omega,~,hyperon]=spin('Omega-');
result=test_true(result,'hyperon angular momentum',...
                 hyperon.spin==1.5&&abs(omega+6.44974765394e7)<1,...
                 'PDG Omega spin is three-halves, not one-half');
result=test_true(result,'hyperon quadrupole unknown',...
                 isnan(hyperon.quad_moment)&&isnan(data('antiOmega-',:).quad_moment),...
                 'spin-three-halves permits a quadrupole moment; unmeasured values are not zero');
result=test_true(result,'particle g-factor convention',...
                 all(isnan(data.gn(data.kind=="particle"))),...
                 'nuclear g factors are inapplicable to particle rows');

% Distinguish unrecognised keys from known but unmeasured signed moments
bad_keys={'not_an_isotope','139Ce','231Am'};
expected={'unknown isotope','unavailable','unavailable'};
for n=1:numel(bad_keys)
    caught=false();
    try
        spin(bad_keys{n});
    catch fault
        caught=contains(fault.message,expected{n});
    end
    result=test_true(result,['diagnostic ' bad_keys{n}],caught,...
                     'missing physical input is a qualified row, not an invented value');
end
result=test_true(result,'unsigned moment',isnan(data('139Ce',:).gamma),...
                 'an experimentally unsigned moment is not given a positive default');

% Exclude subsecond nuclear states without guessing unresolved lifetimes
is_nucleus=data.kind=="nuclide";
subsecond=is_nucleus&data.half_life<1;
result=test_true(result,'isotope lifetime cutoff',...
                 all(data.half_life_kind(subsecond)=="lower limit"),...
                 'a lower bound below one second does not establish a subsecond lifetime');
result=test_true(result,'one-second boundary',data('128Pm',:).half_life==1,...
                 'exactly one second is not shorter than one second');
for name={'6He','229Th_m'}
    caught=false();
    try
        spin(name{1});
    catch fault
        caught=strcmp(fault.identifier,'spin:unknown_isotope');
    end
    result=test_true(result,['excluded ' name{1}],caught,...
                     'the adopted lifetime determines inclusion, not alternative environments');
end

% Preserve synthetic modes and isolate returned metadata from the cache
[gamma,mult,abstract]=spin('G');
result=test_true(result,'ghost metadata',gamma==0&&mult==1&&height(abstract)==0,...
                 'ghosts are not physical isotope records');
for name={'C3','V4','T5','E4'}
    [gamma,mult,abstract]=spin(name{1});
    result=test_true(result,['mode ' name{1}],...
                     mult==str2double(name{1}(2:end))&&height(abstract)==0&&...
                     (gamma==0||gamma==elec),...
                     'synthetic levels preserve the existing numerical contract');
end
[gamma,~,meta]=spin('13C'); meta.gamma=0;
[again,~,meta]=spin('13C');
result=test_true(result,'metadata value semantics',again==gamma&&meta.gamma==gamma,...
                 'editing a returned table must not mutate the persistent database');

% Compare a spin-zero-plus-proton Hilbert system with the proton alone
sys.magnet=1; sys.isotopes={'12C','1H'};
inter.zeeman.scalar={0,1}; bas.formalism='zeeman-hilb'; bas.approximation='none';
spin_zero=test_spin_system(sys,inter,bas);
sys.isotopes={'1H'}; inter.zeeman.scalar={1};
spin_half=test_spin_system(sys,inter,bas);
result=test_close(result,'spin-zero integration',...
                  hamiltonian(assume(spin_zero,'nmr')),...
                  hamiltonian(assume(spin_half,'nmr')),1e-12,1e-12,...
                  'a spin-zero factor cannot change the magnetic Hamiltonian');

end


