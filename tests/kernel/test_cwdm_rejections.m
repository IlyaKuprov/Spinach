% Tests rejection of cross-substance physical specifications. Syntax:
%
%                       result=test_cwdm_rejections()
%
% Outputs:
%
%     result - T3 checks of coupling, product, and symmetry rejections
%
% Direct-sum substances cannot carry pair interactions or product operators
% between blocks. Symmetry inputs use local columns and reject an index
% outside their substance. Error identifiers and messages are asserted.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_rejections()

% Announce the specification boundary
fprintf('TESTING: CWDM cross-substance rejections (T3)\n');
result=new_test_result('kernel/cwdm_rejections','CWDM rejections',...
                      'Cross-substance physics is an error, never a zero.');

% Construct independent two-spin substances
sys.magnet=1; sys.isotopes={'1H','1H','1H','1H'};
inter.chem.parts={1:2,3:4}; inter.chem.concs=[1 0];
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
s=test_spin_system(sys,inter,bas);

% Reject a scalar coupling between different substances
bad=inter; bad.coupling.scalar=cell(4); bad.coupling.scalar{2,3}=10;
rejected=false;
try
    test_spin_system(sys,bad,bas);
catch err
    rejected=strcmp(err.identifier,'Spinach:create:crossSubstanceCoupling')&&...
             contains(err.message,'couplings detected between spins in different chemical species');
end
result=test_true(result,'cross coupling',rejected,...
                 'create raises the named cross-substance coupling error');

% Reject products even when one explicitly requested factor is identity
for labels={{'Lz','Lz'},{'E','Lz'}}
    for constructor={@operator,@state}
        rejected=false;
        try
            constructor{1}(s,labels{1},{1,3});
        catch err
            rejected=strcmp(err.identifier,'Spinach:which_subst:crossSubstance')&&...
                     contains(err.message,'spin list crosses chemical boundaries');
        end
        result=test_true(result,['product ' func2str(constructor{1}) ' ' labels{1}{1}],...
                         rejected,'cross-substance product specifications raise a named error');
    end
end

% Reject a symmetry label that would address another substance globally
bad=bas; bad.sym_group={{'S2'},{}}; bad.sym_spins={{[1 3]},{}};
rejected=false;
try
    basis(s,bad);
catch err
    rejected=contains(err.message,'incorrect spin labels in bas.sym_spins');
end
result=test_true(result,'cross symmetry',rejected,...
                 'symmetry indices must be local to the declared substance');

% Segmented Zeeman filters reject before any tensor-product construction
sys.isotopes={'1H','1H'}; inter.chem.parts={1,2};
for formalism={'zeeman-hilb','zeeman-liouv'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas);
    rho=ones(s.bas.offsets(end),1);
    if strcmp(formalism{1},'zeeman-hilb'), rho=diag(rho); end
    for selector={'correlation','decouple','homospoil'}
        rejected=false;
        try
            switch selector{1}
                case 'correlation'
                    correlation(s,rho,0,1);
                case 'decouple'
                    [~,rho]=decouple(s,[],rho,1);
                case 'homospoil'
                    homospoil(s,rho,'destroy');
            end
        catch err
            rejected=strcmp(err.identifier,['Spinach:' selector{1} ':segmentedZeeman']);
        end
        result=test_true(result,[selector{1} ' ' formalism{1}],rejected,...
                         'unsupported direct-sum Zeeman filtering raises the named error');
    end
    bad=bas; bad.sym_group={{'S2'},{}}; bad.sym_spins={{1},{}};
    rejected=false;
    try
        basis(s,bad);
    catch err
        rejected=strcmp(err.identifier,'Spinach:basis:segmentedZeeman');
    end
    result=test_true(result,['Zeeman symmetry rejection ' formalism{1}],rejected,...
                     'only genuinely segmented Zeeman symmetry is rejected');
end

end
