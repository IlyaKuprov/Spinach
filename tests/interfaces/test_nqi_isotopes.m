% Test: interfaces/nqi_isotopes | Quadrupolar and isotope label import
% Tests Gaussian isotope labels and quadrupolar tensor import. Syntax:
%
%                       result=test_nqi_isotopes()
%
% Outputs:
%
%     result - regression checks for default and explicit isotope labels,
%              the spin-1 tensor convention, isotope substitution, and
%              refusal of ambiguous parser output
%
% talos@spindynamics.org

function result=test_nqi_isotopes()

% Start the test record
result=new_test_result('interfaces/nqi_isotopes',...
                       'Quadrupolar and isotope label import',...
                       'import uses the stated isotope or the magnetic isotope assumed by Gaussian');

% Unlabelled benzoquinone and ammonium carry the Gaussian magnetic isotopes
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
log_dir=fullfile(root_dir,'examples','standard_systems');
log_path=fullfile(log_dir,'benzoquinone.log');
props=gparse(log_path);
expected={'13C','13C','13C','1H','1H','17O','13C','13C','13C','1H','17O','1H'};
[sys,inter_ref]=g2spinach(props,{'C','H','O'},[0 0 0]);
result=test_true(result,'unlabelled benzoquinone',isequal(props.isotopes(:)',expected)&&isequal(sys.isotopes,expected),...
                 'Gaussian assumes 13C, 1H, and 17O for these atoms in coordinate order');
props=gparse(fullfile(log_dir,'shielding_reference','ammonium.log'));
[sys,~]=g2spinach(props,{'N'},0);
result=test_true(result,'unlabelled ammonium',isequal(sys.isotopes,{'14N'}),...
                 'Gaussian assumes 14N for nitrogen');

% Label only the fourth atom in a temporary copy of the benzoquinone input specification
log_text=regexprep(fileread(log_path),'H {21}0\.8484','H(Iso=2)              0.8484','once');
log_copy=[tempname '.log'];
file_id=fopen(log_copy,'w');
fwrite(file_id,log_text);
fclose(file_id);
cleaner=onCleanup(@()delete(log_copy));
props=gparse(log_copy);
expected{4}='2H';
result=test_true(result,'single stated isotope',isequal(find(~isnan(props.mass_numbers)),4)&&...
                 (props.mass_numbers(4)==2)&&isequal(props.isotopes(:)',expected),...
                 'only the labelled hydrogen changes from 1H to 2H');

% The copy keeps J-couplings computed for 1H, which contradict the deuterium label
result=test_true(result,'printed J mismatch',local_refused(props,{'C','H','O'},[0 0 0],'props.j_couplings'),...
                 'J-couplings that K-couplings and the stated isotopes do not reproduce are refused');
[sys,inter]=g2spinach(rmfield(props,'j_couplings'),{'C','H','O'},[0 0 0]);
result=test_true(result,'J-couplings removed',isequal(sys.isotopes,expected),...
                 'without printed J-couplings the stated isotopes are accepted');
result=test_close(result,'J rebuilt for deuterium',inter.coupling.scalar{4,1},...
                  (spin('2H')/spin('1H'))*inter_ref.coupling.scalar{4,1},0,1e-12,...
                  'the K-coupling gives the deuterium J-coupling through the gamma ratio');

% Gaussian prints this 14N quadrupole tensor in MHz in the principal axis frame
props=gparse(fullfile(log_dir,'methylaziridine_a.log'));
nitrogen=strcmp(props.symbols,'N');
printed=[-1.2932  0.6251  1.8700;...
          0.6251  1.7170 -2.3127;...
          1.8700 -2.3127 -0.4238];
[sys,inter]=g2spinach(props,{'N'},0);
result=test_true(result,'quadrupole table isotopes',isequal(sys.isotopes,{'14N'})&&...
                 all(strcmp(props.isotopes(strcmp(props.symbols,'C')),'13C')),...
                 'the table names N(14) and C(12); spin-zero carbon-12 leaves the magnetic carbon-13');
result=test_close(result,'standard frame tensor',inter.coupling.matrix{1,1},props.nqi{nitrogen}/2,0,0,...
                  '2I(2I-1) is two for nitrogen-14');
result=test_close(result,'printed principal values',sort(eig(inter.coupling.matrix{1,1})),...
                  1e6*sort(eig(printed))/2,0,1e-12,'rotation to standard orientation preserves the three principal values');

% Parser output that cannot be interpreted without a guess is refused
ambiguous=props; ambiguous.isotopes{nitrogen}='15N';
result=test_true(result,'quadrupole on spin-half',local_refused(ambiguous,{'N'},0,'props.nqi'),...
                 'a quadrupolar tensor on an isotope without a quadrupole moment is refused');
ambiguous=props; ambiguous.isotopes(strcmp(props.symbols,'C'))={'12C'};
result=test_true(result,'non-magnetic isotope',local_refused(ambiguous,{'C'},0,'non-magnetic'),...
                 'a spin-zero isotope cannot be imported');
result=test_true(result,'missing isotopes',local_refused(rmfield(props,'isotopes'),{'N'},0,'props.isotopes'),...
                 'parser output without isotopes is refused');
orca_props=oparse(fullfile(root_dir,'examples','visualisation','porphyrine.out'));
result=test_true(result,'ORCA hydrogen tensors',local_refused(orca_props,{'H'},0,'props.nqi'),...
                 'ORCA labels hydrogens 1H and prints deuterium quadrupolar tensors for them');

% Nitrogen-15 has no quadrupolar interaction
[sys,inter]=isoswap(sys,inter,1,'15N');
result=test_true(result,'quadrupole removal',isequal(sys.isotopes,{'15N'})&&isempty(inter.coupling.matrix{1,1}),...
                 'substituting spin-half nitrogen removes the quadratic self-coupling');

end

% Reports whether g2spinach refuses the input with an error naming the text
function refused=local_refused(props,elements,references,text)
refused=false;
try
    g2spinach(props,elements,references);
catch exception
    refused=contains(exception.message,text);
end
end

% Trust, but verify.
%
% Russian proverb


