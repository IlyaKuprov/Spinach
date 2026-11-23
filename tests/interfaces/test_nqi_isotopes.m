% Test: interfaces/nqi_isotopes | Quadrupolar and isotope label import
%
% The electric field gradient at a nucleus belongs to the molecule; the
% quadrupolar coupling is its product with the quadrupole moment of the
% isotope. A coupling imported for one isotope from a calculation that
% used another must therefore scale by the ratio of their moments. The
% chlorine ratio is checked against two Gaussian 16 calculations of LiCl
% (B3LYP/6-31G(d), output=pickett), which print -2.7154 MHz for Cl(35)
% and -2.1401 MHz for Cl(37). Isotopes named in a molecule specification
% must be reported by the parser without #p output.
%
% talos@spindynamics.org

function result=test_nqi_isotopes()

% Start the test record
result=new_test_result('interfaces/nqi_isotopes',...
                       'Quadrupolar and isotope label import',...
                       'quadrupolar couplings follow the requested isotope');

% Parse a Gaussian log that carries a nitrogen-14 quadrupole tensor
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
log_dir=fullfile(root_dir,'examples','standard_systems');
props=gparse(fullfile(log_dir,'methylaziridine_a.log'));
nitrogen=find(strcmp(props.symbols,'N'));
proton=find(strcmp(props.symbols,'H'),1);
options.no_xyz=1;

% Gradients exist only where the calculation isotope has a quadrupole moment
result=test_true(result,'gradient provenance',...
                 (~isempty(props.efg{nitrogen}))&&isempty(props.efg{proton}),...
                 'nitrogen-14 has a quadrupole moment, the proton has none');

% Same-isotope import returns the tensor that Gaussian prints
[~,inter]=g2spinach(props,{{'N','14N'}},0,options);
result=test_close(result,'same isotope',inter.coupling.matrix{1,1},...
                  props.nqi{nitrogen}/2,0,1e-12,...
                  'spin 1: the coupling matrix is the printed tensor over 2I(2I-1)');

% A quadrupolar isotope without a gradient is reported and not imported
[report,~,inter]=evalc('g2spinach(props,{{''H'',''2H''}},0,options);');
result=test_true(result,'missing gradient',...
                 contains(report,'no electric field gradient')&&...
                 all(cellfun(@isempty,inter.coupling.matrix),'all'),...
                 'the calculation used protons, so deuterium has no data');

% One chlorine gradient imported as both stable isotopes
chloride.symbols={'Cl'};
chloride.efg={diag([1 1 -2])};
[~,light]=g2spinach(chloride,{{'Cl','35Cl'}},0,options);
[~,heavy]=g2spinach(chloride,{{'Cl','37Cl'}},0,options);
result=test_close(result,'chlorine isotope ratio',...
                  light.coupling.matrix{1,1}(3,3)/heavy.coupling.matrix{1,1}(3,3),...
                  2.7154/2.1401,0,1e-4,...
                  'both isotopes have spin 3/2; Gaussian prints five significant digits');

% Label one hydrogen of a shipped log as deuterium and parse the copy
log_text=fileread(fullfile(log_dir,'benzoquinone.log'));
log_text=regexprep(log_text,'H {21}0\.8484','H(Iso=2)              0.8484','once');
log_copy=[tempname '.log'];
file_id=fopen(log_copy,'w');
fwrite(file_id,log_text);
fclose(file_id);
props=gparse(log_copy);
delete(log_copy);
stated=find(~isnan(props.mass_numbers));
result=test_true(result,'isotope label',...
                 (numel(props.mass_numbers)==numel(props.symbols))&&...
                 isscalar(stated)&&strcmp(props.symbols{stated},'H')&&...
                 (props.mass_numbers(stated)==2),...
                 'one labelled hydrogen, no statement about the other atoms');

end

% Trust, but verify.
%
% Russian proverb


