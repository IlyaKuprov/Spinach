% Intermolecular exchange as an explicit reaction-record expansion. Syntax:
%
%                      result=test_cwdm_flux()
%
% Outputs:
%
%    result - magnetisation transfer, correlation loss, and population checks
%
% A two-spin molecule exchanges its first spin with a one-spin pool.
% A+B to A+B keeps both species concentrations fixed. Additive arrival
% carries either reactant's internal orders but discards cross-reactant
% products, reproducing magnetisation exchange and loss of correlations
% involving the departing spin. No cross-molecule correlations are stored.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_flux()

% Compile a spin replacement event between a molecule and its pool
fprintf('TESTING: CWDM intermolecular exchange expansion\n');
result=new_test_result('kernel/cwdm_flux','CWDM intermolecular flux',...
                      'Replacement reactions transfer magnetisation and destroy departing-spin correlations.');
sys.magnet=1; sys.isotopes={'1H','1H','1H'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.chem.parts={1:2,3}; inter.chem.concs=[1 1];
inter.chem.reactions={struct('reactants',[1 2],'products',[1 2],...
                            'matching',[1 3;2 2;3 1],'rate',2)};
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
s=basis(create(sys,inter),bas); generator=kinetics(s);
eta=unit_state(s); K=generator(0,eta);

% Assert transfer and conservation of single-spin magnetisation
molecule=coil_state(s,'Lz',1,'exact'); partner=coil_state(s,'Lz',2,'exact');
pool=coil_state(s,'Lz',3,'exact'); correlation=coil_state(s,{'Lz','Lz'},{1,2},'exact');
result=test_close(result,'molecule to pool',K*molecule,2*(pool-molecule),1e-12,0,...
                  'the departing spin enters the pool at the reaction rate');
result=test_close(result,'pool to molecule',K*pool,2*(molecule-pool),1e-12,0,...
                  'the arriving pool spin replaces the molecular spin');
result=test_close(result,'retained spin',K*partner,zeros(size(partner)),1e-12,0,...
                  'orders not involving the exchanged spin are preserved');
result=test_close(result,'departing-spin correlation',K*correlation,-2*correlation,1e-12,0,...
                  'tracing the departing spin destroys its intramolecular correlation');

% Concentrations remain fixed independently of their initial values
rng(13);
for n=1:8
    trial=randn(size(eta)); trial(s.bas.offsets(1:end-1)+1)=rand(2,1);
    result=test_close(result,['population conservation ' num2str(n)],...
                      chem_concs(s,generator(0,trial)*trial),[0 0],1e-12,0,...
                      'reactant and product stoichiometries are identical');
end

% Fixed concentrations permit freezing this additive generator
eta=eta+0.2*molecule+0.3*pool+0.1*correlation;
reference=expm(full(K)*0.7)*eta; current=eta; dt=0.01;
for n=1:70
    current=step(s,{@(t,y)1i*generator(t,y),(n-1)*dt,'RKMK4'},current,dt);
end
result=test_close(result,'frozen additive flux',current,reference,1e-11,0,...
                  'with invariant concentrations the state-dependent additive generator stays constant');
result=test_close(result,'total exchanged magnetisation',...
                  (molecule+pool)'*current,(molecule+pool)'*eta,1e-12,0,...
                  'replacement moves, but does not destroy, the exchanged single-spin magnetisation');

end


