% Tests the explicit chemical reaction input contract. Syntax:
%
%                       result=test_cwdm_records()
%
% Outputs:
%
%    result - reaction-record validation and merging test results
%
% Every new create rejection is checked by identifier and message. Valid
% records include loss, spin-free products, time-dependent rates, repeated
% substances, and named or user-specified spin selectors.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_records()

% Announce the input contract
fprintf('TESTING: CWDM reaction record inputs\n');
result=new_test_result('kernel/cwdm_records','Reaction records',...
                      'Explicit chemistry input validation and index offsets.');

% Construct independently labelled reactant and product spins
sys.magnet=1; sys.isotopes={'E','E','1H','1H','13C'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.chem.parts={1:3,4,5,[]}; inter.chem.concs=[1 0 0 0];
reaction=struct('reactants',1,'products',2,'matching',[3 4],'rate',2);
inter.chem.reactions={reaction};
s=create(sys,inter);
result=test_true(result,'default closure',strcmp(s.chem.reactions{1}.closure,'additive'),...
                 'the declared default closure is additive');
summary_chemistry(s);

% Preserve column-oriented spin membership while formatting its summary
column_inter.chem.parts={(1:5)'};
column_system=create(sys,column_inter); column_system.sys.output=1;
text=evalc('summary_chemistry(column_system);');
result=test_true(result,'column membership summary',...
                 isequal(column_system.chem.parts{1},(1:5)')&&...
                 contains(text,'chemical subsystem 1: spins [1  2  3  4  5]'),...
                 'reporting formats a row without changing the accepted column-oriented input');

% Assert every retired input even when its value is empty
for field={'rates','flux_rate','flux_type','rp_theory','rp_rates','rp_electrons'}
    bad=inter; bad.chem.(field{1})=[]; rejected=false;
    try
        create(sys,bad);
    catch err
        rejected=strcmp(err.identifier,'Spinach:create:retiredChemistry')&&...
                 contains(err.message,['inter.chem.' field{1} ' is retired'])&&...
                 contains(err.message,'inter.chem.reactions records');
    end
    result=test_true(result,['retired ' field{1}],rejected,'retired fields name the reaction-record replacement');
end

% Enumerate malformed records and the precise rejection contracts
cases={...
    'reactants',[0 1],'reactionSubstances','row vectors of valid substance indices';...
    'reactants',[],'reactionReactants','at least one reactant';...
    'matching',[3 4 5],'reactionMatching','two-column matrix';...
    'matching',[4 3],'reactionMembership','declared reactants and products';...
    'matching',[3 4;3 4],'reactionDuplicate','must not be matched twice';...
    'products',3,'reactionMembership','declared reactants and products';...
    'rate',-1,'reactionRate','finite non-negative scalar';...
    'closure','other','reactionClosure','additive or product';...
    'selector',{'singlet'},'reactionSelector','two-element cell';...
    'selector',{'singlet',[1 3]},'reactionElectrons','two distinct electrons';...
    'selector',{ones(2),ones(3)},'reactionSelectorPair','equal size';...
    'extra',1,'reactionField','unrecognised reaction record field'};
for n=1:size(cases,1)
    bad=inter; bad.chem.reactions{1}.(cases{n,1})=cases{n,2}; rejected=false;
    try
        create(sys,bad);
    catch err
        rejected=strcmp(err.identifier,['Spinach:create:' cases{n,3}])&&contains(err.message,cases{n,4});
    end
    result=test_true(result,['record ' num2str(n)],rejected,'identifier and explanatory message match');
end

% Exercise container, required-field, partition, and isotope errors
for n=1:5
    bad=inter;
    switch n
        case 1
            bad.chem.reactions=1; id='reactionRecords'; text='cell vector';
        case 2
            bad.chem.reactions={rmfield(reaction,'rate')}; id='reactionRecord'; text='rate fields';
        case 3
            bad.chem=rmfield(bad.chem,'parts'); id='reactionParts'; text='explicit inter.chem.parts';
        case 4
            bad.chem.extra=1; id='chemistryField'; text='only parts, concs, and reactions';
        case 5
            bad.chem.reactions{1}.products=3; bad.chem.reactions{1}.matching=[3 5];
            id='reactionIsotopes'; text='identical isotopes';
    end
    rejected=false;
    try
        create(sys,bad);
    catch err
        rejected=strcmp(err.identifier,['Spinach:create:' id])&&contains(err.message,text);
    end
    result=test_true(result,id,rejected,'identifier and explanatory message match');
end

% Accept loss, sinks, multiplicities, time dependence, and selectors
valid={reaction,reaction,reaction,reaction,reaction};
valid{1}.products=[]; valid{1}.matching=zeros(0,2);
valid{2}.products=4; valid{2}.matching=zeros(0,2); valid{2}.rate=@(t)2+t;
valid{3}.reactants=[1 1]; valid{3}.products=[2 2]; valid{3}.closure='product';
valid{4}.selector={'singlet',[1 2]}; valid{5}.selector={speye(64),speye(64)};
inter.chem.reactions=valid; s=create(sys,inter);
result=test_true(result,'valid records',numel(s.chem.reactions)==5,'valid record variants are retained');

% Preserve the default single substance when chemistry groups are empty
empty.chem=struct();
[empty_sys,empty_inter]=merge_inp({sys,sys},{empty,empty});
empty_system=create(empty_sys,empty_inter);
result=test_true(result,'empty chemistry merge',isempty(fieldnames(empty_inter.chem))&&...
                 isequal(empty_system.chem.parts,{1:10})&&isempty(empty_system.chem.reactions),...
                 'merging absent partitions must not invent a reaction field requiring explicit parts');

% Merge independent records with substance and spin offsets
[merged_sys,merged]=merge_inp({sys,sys},{inter,inter});
s=create(merged_sys,merged);
result=test_true(result,'merged matching',isequal(s.chem.reactions{6}.reactants,5)&&...
                 isequal(s.chem.reactions{9}.matching,[8 9])&&...
                 isequal(s.chem.reactions{9}.selector{2},[6 7])&&...
                 isequal(s.chem.reactions{10}.selector{1},speye(64)),...
                 'substance and named spin indices shift; local selector matrices do not');

end


