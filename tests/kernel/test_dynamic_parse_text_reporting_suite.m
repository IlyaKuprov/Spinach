% Tests deterministic parsing, text, and safe reporting utilities. Syntax:
%
%                    result=test_dynamic_parse_text_reporting_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks operator-specification parsing, isotope predicates,
% label lookup, silent reporting calls, and polyadic text diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_parse_text_reporting_suite()

% Announce the test target
fprintf('TESTING: Parsing, text, and reporting utilities\n');

% State the utility target of the test
result=new_test_result('kernel/dynamic_parse_text_reporting_suite',...
                       'Parsing, text, and reporting utilities',...
                       'deterministic parsing and silent reporting helpers must preserve documented text and metadata semantics.');

% Build a small spin-system descriptor for parsing helpers
spin_system=local_parse_spin_system();

% Check isotope and spin-label parsing into single-spin opspecs
[opspecs,coeffs]=human2opspec(spin_system,'Lz','nuclei');
result=test_true(result,'human2opspec nuclei Lz',...
                 isequal(opspecs,{[2 0 0];[0 0 2]})&&isequal(coeffs,[1;1]),...
                 'nuclei must select non-electron spins and Lz must map to IST index two');

% Check product-operator parsing and Lx expansion coefficients
[prod_ops,prod_coeffs]=human2opspec(spin_system,{'Lx','Lz'},{1,3});
result=test_true(result,'human2opspec product opspecs',...
                 isequal(prod_ops,{[1 0 2];[3 0 2]}),...
                 'Lx on spin one expands into L+ and L- terms while Lz remains on spin three');
result=test_close(result,'human2opspec Lx coefficients',prod_coeffs,[-sqrt(2);sqrt(2)]/2,1e-15,1e-15,...
                  'Spinach spherical tensor convention gives Lx=(L+ + L-)/2');

% Check label lookup against a unique label list
sys.isotopes=spin_system.comp.isotopes;
sys.labels=spin_system.comp.labels;
result=test_close(result,'idxof label lookup',idxof(sys,'carbon'),3,1e-15,1e-15,...
                  'idxof must return the one-based spin index matching the requested label');

% Check electron and nucleus isotope predicates
result=test_true(result,'isnucleus proton',isnucleus('1H'),...
                 'ordinary isotope specifications must be classified as nuclei');
result=test_true(result,'iselectron electron',iselectron('E'),...
                 'the electron particle specification must be classified as an electron');
result=test_true(result,'isnucleus electron false',~isnucleus('E'),...
                 'electron specifications must not be classified as nuclei');

% Check that hush-mode report is side-effect-free on the console
report_text=evalc('report(spin_system,''hidden message'');');
result=test_true(result,'report hush output',isempty(report_text),...
                 'report must produce no console text when spin_system.sys.output is hush');

% Check that banner delegates safely through hush-mode report
banner_text=evalc('banner(spin_system,''basis_banner'');');
result=test_true(result,'banner hush output',isempty(banner_text),...
                 'banner must complete silently when reporting is hushed');

% Check that summary can traverse coordinate metadata without printing
summary_text=evalc('summary_coordinates(spin_system,''coordinate summary'');');
result=test_true(result,'summary hush output',isempty(summary_text),...
                 'summary_coordinates must respect hush-mode reporting while reading coordinate metadata');

% Check polyadic text diagnostics on a small unopened Kronecker product
info_text=evalc('polinfo(polyadic({{speye(2),sparse(1)}}));');
result=test_true(result,'polinfo shape output',contains(info_text,'polyadic [2x2]')&&...
                 contains(info_text,'matrix 1 [2x2]')&&contains(info_text,'matrix 2 [1x1]'),...
                 'polinfo must describe the unopened polyadic size and its matrix cores');

end


function spin_system=local_parse_spin_system()

% Create quiet system output settings
spin_system.sys.output='hush';
spin_system.sys.disable={};

% Define three simple spin particles
spin_system.comp.nspins=3;
spin_system.comp.isotopes={'1H','E','13C'};
spin_system.comp.labels={'proton','electron','carbon'};
spin_system.comp.types={'S','S','S'};
spin_system.comp.mults=[2 2 2];

% Attach coordinates for summary traversal
spin_system.inter.coordinates={[0 0 0],[1 0 0],[0 1 0]};

end


