% Tests direct-sum Zeeman dimensions and storage-only wavefunction limits.
% Syntax:
%
%                      result=test_cwdm_formalisms()
%
% Outputs:
%
%    result - exact dimension, offset, and named capability checks
%
% Unequal Hilbert dimensions and a spin-free block distinguish a direct
% sum from a tensor product. Unsupported wavefunction operations assert
% both the error identifier and the complete formalism-naming message.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_formalisms()

% Announce the formalism contract
result=new_test_result('kernel/cwdm_formalisms','CWDM formalism capabilities',...
                       'Zeeman storage uses local dimensions and explicit capability limits.');

% Reject absent basis metadata before inspecting formalism capabilities
for incomplete={struct(),struct('bas',struct())}
    rejected=false;
    try
        unit_state(incomplete{1});
    catch err
        rejected=strcmp(err.message,'the spin_system object does not contain the required information.');
    end
    result=test_true(result,'unit basis metadata',rejected,...
                     'missing basis or formalism retains the explicit input-validation error');
end

% Include unequal spin-bearing blocks and a spin-free substance
sys.magnet=14.1; sys.isotopes={'1H','1H','1H'};
inter.temperature=298; inter.zeeman.scalar={1,2,3};
inter.chem.parts={1,2:3,[]}; inter.chem.concs=[0.4 0.6 0];
bas.approximation={'none','none','none'};
for formalism={'zeeman-liouv','zeeman-hilb','zeeman-wavef'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas);
    dims=[2;4;1];
    if strcmp(formalism{1},'zeeman-liouv'), dims=dims.^2; end
    result=test_true(result,['dimensions ' formalism{1}],isequal(s.bas.nstates,dims),...
                     'each block has its local Hilbert dimension or its square');
    result=test_true(result,['offsets ' formalism{1}],isequal(s.bas.offsets,[0;cumsum(dims)]),...
                     'direct-sum offsets include the spin-free scalar block');
    result=test_true(result,['descriptors ' formalism{1}],all(cellfun(@isempty,s.bas.basis)),...
                     'Zeeman blocks do not store spherical-tensor descriptors');
    fprintf('CWDM_FORMALISM_DIM %s %s\n',formalism{1},mat2str(dims'));
end

% Store local pure kets without assigning concentrations to their amplitudes
ket=coil_state(s,[0.5 0.5 0.5],[],'exact');
result=test_close(result,'unweighted wavefunction storage',ket,[1;0;1;0;0;0;1],0,0,...
                  'each substance stores its own product ket, including the scalar empty block');
rejected=false;
try
    state(s,[0.5 0.5 0.5]);
catch err
    rejected=strcmp(err.identifier,'Spinach:state:segmentedZeeman')&&...
             strcmp(err.message,'concentration-weighted states are not supported in segmented zeeman-wavef formalism.');
end
result=test_true(result,'wavefunction concentration weighting',rejected,...
                 'the weighted wrapper does not reinterpret pure-state amplitudes as concentrations');

% Check thermal and concentration-weighted wavefunction requests explicitly
calls={@()equilibrium(s),@()unit_state(s),...
       @()thermalize(s,sparse(7,7),[],[],[], 'IME'),...
       @()thermalize(s,sparse(7,7),speye(7),298,[], 'dibari'),...
       @()chem_concs(s,ket)};
ids={'Spinach:equilibrium:wavefunction','Spinach:unit_state:wavefunction',...
     'Spinach:thermalize:wavefunction','Spinach:thermalize:wavefunction',...
     'Spinach:chem_concs:formalism'};
messages={'thermal equilibrium is not supported in zeeman-wavef formalism.',...
          'concentration-weighted unit states are not supported in zeeman-wavef formalism.',...
          'thermalisation is not supported in zeeman-wavef formalism.',...
          'thermalisation is not supported in zeeman-wavef formalism.',...
          'concentrations are not defined in zeeman-wavef formalism.'};
for n=1:numel(calls)
    rejected=false;
    try
        calls{n}();
    catch err
        rejected=strcmp(err.identifier,ids{n})&&strcmp(err.message,messages{n});
    end
    result=test_true(result,ids{n},rejected,'the message names the unsupported wavefunction capability');
end

% Reject first-order and mass-action chemistry at basis and kinetics entry
for reactants={1,[1 2]}
    reaction=struct('reactants',reactants{1},'products',3,...
                    'matching',zeros(0,2),'rate',2,'closure','additive');
    s.chem.reactions={reaction};
    calls={@()basis(s,bas),@()kinetics(s)};
    ids={'Spinach:basis:wavefunctionChemistry','Spinach:kinetics:wavefunction'};
    for n=1:numel(calls)
        rejected=false;
        try
            calls{n}();
        catch err
            rejected=strcmp(err.identifier,ids{n})&&...
                     strcmp(err.message,'chemical reactions are not supported in zeeman-wavef formalism.');
        end
        result=test_true(result,[ids{n} mat2str(reactants{1})],rejected,...
                         'both first-order and mass-action records are rejected by name');
    end
end
fprintf('CWDM_FORMALISM_CAPABILITIES failures=%d\n',numel(result.failures));
end


