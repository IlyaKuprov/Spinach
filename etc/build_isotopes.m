% Build the typed runtime isotope table from the audited literature TSV.
% Syntax:
%
%                  isotopes=build_isotopes(source_file)
%
% Parameters:
%
%    source_file  - character row vector, path to the isotope TSV;
%                   the MAT payload is saved alongside that file
%
% Outputs:
%
%    isotopes     - physical-data table with explicit units, row names,
%                   source qualifiers, and derived nuclear g factors
%
% This is an offline data-preparation operation, never a lookup fallback.
% Source DOI attribution, completeness, and adoption choices must already
% have been audited. Nuclear moments are in nuclear magnetons; quadrupole
% moments are in barns. CODATA2022 mu_N and exact hbar define the conversion
% to angular magnetogyric ratio. Particle gamma values are supplied by their
% independently cited source, not calculated as nuclear g factors.
%
% ilya.kuprov@weizmann.ac.il

function isotopes=build_isotopes(source_file)

% Check consistency
grumble(source_file);

% Import the reviewed schema without text type inference
names={'isotope','kind','state','coverage_status','z','a','element',...
       'element_name','spin','spin_status','mag_moment','gn','gamma',...
       'abundance','abund_low','abund_high','abund_basis','quad_moment',...
       'radioactive','half_life','half_life_kind','source_dois','notes'};
types={'string','string','string','string','double','double','string',...
       'string','double','string','double','double','double','double',...
       'double','double','string','double','logical','double','string',...
       'string','string'};
opts=delimitedTextImportOptions('NumVariables',numel(names));
opts.Delimiter='\t'; opts.VariableNames=names; opts.VariableTypes=types;
opts.DataLines=[2 Inf]; opts.ExtraColumnsRule='error';
isotopes=readtable(source_file,opts);

% Derive nuclear g factors and gamma with one explicit constants convention
mu_n=5.0507837393e-27; hbar=6.62607015e-34/(2*pi);
is_nucleus=isotopes.kind=="nuclide";
nonzero=is_nucleus&isotopes.spin>0;
isotopes.gn(nonzero)=isotopes.mag_moment(nonzero)./isotopes.spin(nonzero);
derived=nonzero&~isnan(isotopes.mag_moment)&isnan(isotopes.gamma);
isotopes.gamma(derived)=isotopes.gn(derived)*mu_n/hbar;
zero_spin=is_nucleus&isotopes.spin==0;
isotopes.mag_moment(zero_spin)=0; isotopes.gamma(zero_spin)=0;
isotopes.quad_moment(is_nucleus&isotopes.spin<1)=0;
isotopes.source_dois(derived)=isotopes.source_dois(derived)+...
    "; gamma,constants: 10.1103/RevModPhys.97.025002";

% Enforce physically meaningful source-independent invariants
assert(numel(unique(isotopes.isotope))==height(isotopes),'duplicate isotope labels.');
assert(all(isotopes.z(is_nucleus)>=1&isotopes.z(is_nucleus)==fix(isotopes.z(is_nucleus))),...
       'nuclear atomic numbers must be positive integers.');
assert(all(isotopes.a(is_nucleus)>=isotopes.z(is_nucleus)&...
           isotopes.a(is_nucleus)==fix(isotopes.a(is_nucleus))),...
       'nuclear mass numbers must be integers no smaller than atomic numbers.');
known=~isnan(isotopes.spin);
assert(all(isotopes.spin(known)>=0&2*isotopes.spin(known)==fix(2*isotopes.spin(known))),...
       'known spins must be nonnegative integer or half-integer.');
for name={'abundance','abund_low','abund_high'}
    values=isotopes.(name{1}); known=~isnan(values);
    assert(all(values(known)>=0&values(known)<=1),'abundance fractions are outside [0,1].');
end
known=~isnan(isotopes.abund_low)&~isnan(isotopes.abund_high);
assert(all(isotopes.abund_low(known)<=isotopes.abund_high(known)),...
       'abundance bounds are reversed.');
assert(all(isnan(isotopes.half_life)|isotopes.half_life>0),...
       'known half-lives must be positive.');
assert(all(contains(isotopes.source_dois,"10.")),'every physical row needs academic DOI attribution.');

% Attach physical units and stable schema metadata before saving
isotopes.Properties.RowNames=cellstr(isotopes.isotope);
units=repmat({''},1,numel(names));
units{11}='nuclear magneton'; units{13}='rad/(s*T)';
units{14}='fraction'; units{15}='fraction'; units{16}='fraction';
units{18}='barn'; units{20}='s';
isotopes.Properties.VariableUnits=units;
isotopes.Properties.VariableDescriptions={...
    'Canonical exact lookup label','Nuclide or particle','Nuclear state or free particle',...
    'Inclusion qualification','Atomic number','Mass number','Element symbol','Element name',...
    'Angular momentum quantum number','Spin assignment qualification',...
    'Signed magnetic dipole moment','Nuclear moment divided by nonzero spin',...
    'Signed angular magnetogyric ratio','Representative terrestrial abundance',...
    'Published terrestrial abundance lower bound','Published terrestrial abundance upper bound',...
    'Abundance convention and evaluation','Signed spectroscopic quadrupole moment',...
    'Finite decay or unclassified lifetime','Adopted half-life in seconds',...
    'Lifetime value, bound, estimate, stability, or unknown',...
    'Field-labelled academic source DOIs','Source uncertainties, conditions, and adoption qualifications'};
isotopes.Properties.Description='Literature-sourced nuclear and particle spin properties';
isotopes.Properties.UserData=struct('schema',1,'constants','CODATA2022',...
    'half_life_note','Inf denotes classified stable; qualifiers and source conditions remain explicit');

% Save an uncompressed release payload without a runtime regeneration path
[folder,stem]=fileparts(source_file);
save(fullfile(folder,[stem '.mat']),'isotopes','-v7','-nocompression');

end

% Consistency enforcement
function grumble(source_file)
if ~ischar(source_file)||isempty(source_file)||~isrow(source_file)
    error('source_file must be a nonempty character row vector.');
end
end


