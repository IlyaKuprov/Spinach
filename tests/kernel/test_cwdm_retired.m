% Tests retired basis and chemistry fields at consumer boundaries. Syntax:
%
%                       result=test_cwdm_retired()
%
% Outputs:
%
%     result - named retirement error and replacement-message checks
%
% Ordinary MATLAB structs cannot intercept arbitrary external dot reads.
% These tests pass legacy layouts to the basis compiler, descriptor filters,
% basis summary, and kinetics consumer; valid segmented inputs also run.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_retired()

% Announce the retirement contract
fprintf('TESTING: CWDM retired consumer fields\n');
result=new_test_result('kernel/cwdm_retired','Retired consumer fields',...
                      'Named errors point from legacy layouts to replacement fields.');

% Construct a valid local descriptor and physical state
sys.magnet=1; sys.isotopes={'1H'};
inter.chem.concs=1;
bas.formalism='sphten-liouv'; bas.approximation={'none'};
s=test_spin_system(sys,inter,bas);
rho=state(s,'Lz','1H');

% Reject both legacy basis fields at each named boundary
for n=1:2
    bad=s; bad_bas=bas;
    if n==1
        bad.bas.basis=s.bas.basis{1}; bad_bas.basis=s.bas.basis{1};
        id='Spinach:basis:retiredGlobalBasis';
        text={'global bas.basis matrix is retired','bas.basis{n}','bas.offsets'};
    else
        bad.bas.irrep=[]; bad_bas.irrep=[];
        id='Spinach:basis:retiredIrrep';
        text={'bas.irrep is retired','bas.sym_fact(n).irr_projectors','irr_dimensions'};
    end
    calls={@()basis(s,bad_bas),@()coherence(bad,rho,{{'1H',0}}),...
           @()correlation(bad,rho,1,'all'),@()summary_basis(bad),@()kinetics(bad)};
    names={'basis','coherence','correlation','summary_basis','kinetics'};
    for k=1:numel(calls)
        rejected=false;
        try
            calls{k}();
        catch err
            rejected=strcmp(err.identifier,id)&&all(contains(err.message,text));
        end
        result=test_true(result,[names{k} ' field ' num2str(n)],rejected,...
                         'the named error states the retired field and its replacements');
    end
end

% Reject chemistry fields inserted after create even when empty
for field={'rates','flux_rate','flux_type','rp_theory','rp_rates','rp_electrons'}
    bad=s; bad.chem.(field{1})=[]; rejected=false;
    try
        kinetics(bad);
    catch err
        rejected=strcmp(err.identifier,'Spinach:kinetics:retiredChemistry')&&...
                 contains(err.message,['chem.' field{1} ' is retired'])&&...
                 contains(err.message,'chem.reactions records');
    end
    result=test_true(result,['retired ' field{1}],rejected,...
                     'the chemistry consumer names the reaction-record replacement');
end

% Retain supported consumers and the normal compiled structure
summary_basis(s); rebuilt=basis(s,bas);
result=test_true(result,'valid compiler',isequal(rebuilt.bas.basis,s.bas.basis),...
                 'ordinary cell descriptors compile without retirement errors');
result=test_close(result,'valid coherence',coherence(s,rho,{{'1H',0}}),rho,0,0,...
                  'longitudinal magnetisation has zero coherence order');
result=test_close(result,'valid correlation',correlation(s,rho,1,'all'),rho,0,0,...
                  'single-spin magnetisation has correlation order one');
result=test_true(result,'valid kinetics',nnz(kinetics(s))==0,...
                 'a chemistry-free system retains its zero generator');
fprintf('CWDM_RETIRED_MESSAGES_COMPLETE basis=10 chemistry=6\n');

end


