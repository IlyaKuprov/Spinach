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
molecule=coil_state(s,'Lz',1); partner=coil_state(s,'Lz',2);
pool=coil_state(s,'Lz',3); correlation=coil_state(s,{'Lz','Lz'},{1,2});
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

% Original point-to-point flux syntax remains a distinct additive R contribution
result=magnetisation_flux_checks(result);
result=flux_consumer_checks(result);
result=inactive_flux_checks(result);

end

% Check the historical point-to-point model independently of reaction records
function result=magnetisation_flux_checks(result)

% Retain the original one-substance spin space and non-zero diagonal bath rates
sys.magnet=14.1; sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={0,3,5};
inter.chem.flux_rate=[11 500 0;2000 13 70;0 30 17];
inter.chem.flux_type='intermolecular';
bas.formalism='sphten-liouv'; bas.approximation={'none'};
s=test_spin_system(sys,inter,bas);

% Spin removal retains the original flux matrix indexing contract
reduced=kill_spin(s,3);
result=test_close(result,'flux indices after spin removal',...
                  reduced.chem.flux_rate,inter.chem.flux_rate(1:2,1:2),0,0,...
                  'removed spins lose their incoming and outgoing flux entries');

result=test_true(result,'original flux inputs retained',...
                 isequal(s.chem.flux_rate,inter.chem.flux_rate)&&...
                 strcmp(s.chem.flux_type,inter.chem.flux_type),...
                 'point-to-point rates retain their original matrix orientation');
for type={'intermolecular','intramolecular'}
    s.chem.flux_type=type{1};
    result=test_close(result,[type{1} ' full flux generator'],...
                      relaxation(s),original_flux(s),0,0,...
                      'single-spin gains, losses, correlations, and diagonal rates match the original generator');
end
result=test_close(result,'flux not duplicated in kinetics',kinetics(s),...
                  sparse(s.bas.offsets(end),s.bas.offsets(end)),0,0,...
                  'point-to-point flux belongs only to R');

% Local descriptors and non-contiguous spin labels must map to the correct blocks
multi=inter; multi.chem.parts={[1 3],2}; multi.chem.concs=[1 1];
local=bas; local.approximation={'none','none'};
m=test_spin_system(sys,multi,local);
result=test_close(result,'flux across local descriptor blocks',...
                  relaxation(m),original_flux(m),0,0,...
                  'global spin labels are translated to local descriptor columns and offsets');

% Relaxation retention and thermalisation must not act on the flux generator
s.chem.flux_type='intermolecular';
for keep={'diagonal','secular','labframe'}
    for destination={'zero','IME','dibari'}
        work=s; work.rlx.theories={'t1_t2'};
        work.rlx.r1_rates={1,2,3}; work.rlx.r2_rates={4,5,6};
        work.rlx.keep=keep{1}; work.rlx.equilibrium=destination{1};
        work.rlx.temperature=298;
        if strcmp(destination{1},'dibari')
            work.sys.disable=[work.sys.disable {'high_temp'}];
        end
        bare=work; bare.chem.flux_rate=[]; bare.chem.flux_type='';
        result=test_close(result,['flux after ' keep{1} ' ' destination{1}],...
                          relaxation(work)-relaxation(bare),original_flux(work),1e-12,1e-12,...
                          'flux is added after relaxation processing, as the old separate K was');
    end
end
work.tols.rlx_zero=1e5;
result=test_close(result,'flux bypasses relaxation zero threshold',...
                  relaxation(work),original_flux(work),0,0,...
                  'relaxation-specific numerical thresholds must not erase historical flux');

% Original no-flux defaults and paired-field validation remain unchanged
bare=inter; bare=rmfield(bare,'chem');
bare=test_spin_system(sys,bare,bas);
result=test_true(result,'original no-flux defaults',...
                 isempty(bare.chem.flux_rate)&&isempty(bare.chem.flux_type),...
                 'absent flux inputs produce the original empty defaults');
invalid={struct('flux_rate',zeros(3)),struct('flux_type','intermolecular'),...
         struct('flux_rate',zeros(2),'flux_type','intermolecular'),...
         struct('flux_rate',1i*eye(3),'flux_type','intermolecular'),...
         struct('flux_rate',zeros(3),'flux_type','unknown')};
for n=1:numel(invalid)
    bad=inter; bad.chem=invalid{n}; rejected=false;
    try
        create(sys,bad);
    catch
        rejected=true;
    end
    result=test_true(result,['original flux input validation ' num2str(n)],rejected,...
                     'the paired fields, real square spin-sized matrix, and flux type are required');
end

% Unlike tensor ranks cannot be transferred by the original basis matching
sys.isotopes={'1H','2H'}; inter.zeeman.scalar={0,3};
inter.chem.flux_rate=[0 3;0 0];
s=test_spin_system(sys,inter,bas); rejected=false;
try
    relaxation(s);
catch err
    rejected=contains(err.message,'different topology or basis sets');
end
result=test_true(result,'original tensor-label compatibility',rejected,...
                 'different local tensor ranks retain the original basis-mismatch rejection');

end

% Historical row-mask definition, evaluated on the current direct-sum row labels
function F=original_flux(s)

% Embed local descriptors only in this small independent regression reference
b=zeros(s.bas.offsets(end),s.comp.nspins);
for n=1:numel(s.chem.parts)
    rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    b(rows,s.chem.parts{n})=s.bas.basis{n};
end
orders=sum(logical(b),2); F=sparse(size(b,1),size(b,1));
[from,to,rates]=find(s.chem.flux_rate);
for n=1:numel(rates)
    for order=1:2
        mask=(orders==1);
        if order==2, mask=(orders>1); end
        source=mask&(b(:,from(n))~=0);
        destin=mask&(b(:,to(n))~=0);
        if (order==1)||strcmp(s.chem.flux_type,'intramolecular')
            stationary=source&destin;
            source=xor(source,stationary); destin=xor(destin,stationary);
            if ~isequal(b(source,from(n)),b(destin,to(n)))
                error('reference tensor labels do not match.');
            end
            F=F+rates(n)*sparse(find(destin),find(source),ones(nnz(source),1),size(F,1),size(F,2));
        end
        F=F-rates(n)*spdiags(double(source),0,size(F,1),size(F,2));
    end
end

end



% Flux must remain compatible with recursive relaxation, assembly, and propagation
function result=flux_consumer_checks(result)

% Use a non-zero, valid SRSK coupling and flux between two proton sites
sys.magnet=14.1; sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={0,3,5}; inter.coupling.scalar=cell(3);
inter.coupling.scalar{1,3}=3;
inter.relaxation={'t1_t2','SRSK'}; inter.srsk_sources=3;
inter.r1_rates={1,2,1e5}; inter.r2_rates={3,4,1e5};
inter.rlx_keep='labframe'; inter.equilibrium='zero'; inter.temperature=298;
inter.chem.flux_rate=[11 5 0;7 13 0;0 0 0];
inter.chem.flux_type='intermolecular';
bas.formalism='sphten-liouv'; bas.approximation={'none'};
s=test_spin_system(sys,inter,bas);
flux=s; flux.rlx.theories={}; F=relaxation(flux);
for destination={'zero','IME','dibari'}
    work=s; work.rlx.equilibrium=destination{1};
    if strcmp(destination{1},'dibari')
        work.sys.disable=[work.sys.disable {'high_temp'}];
    end
    bare=work; bare.chem.flux_rate=[]; bare.chem.flux_type='';
    bath=relaxation(bare);
    result=test_close(result,['SRSK once-only flux ' destination{1}],...
                      relaxation(work)-bath,F,100*eps(norm(bath,'fro')),0,...
                      'recursive SRSK rates exclude flux until after outer thermalisation');
end

% Reject constituent sizes even when the merged dimensions would match
one.isotopes={'1H'}; two.isotopes={'1H','1H'};
a.chem.flux_rate=[0 1;2 0]; a.chem.flux_type='intermolecular';
b.chem.flux_rate=0; b.chem.flux_type='intermolecular'; rejected=false;
try
    merge_inp({one,two},{a,b});
catch err
    rejected=contains(err.message,'subsystem spin count');
end
result=test_true(result,'compensating flux block sizes rejected',rejected,...
                 'each block must match its own subsystem before concatenation');

% Reuse the historical block-diagonal and common-value merge contracts
part.isotopes={'1H','1H'};
a.chem.flux_rate=[3 -5;7 0]; a.chem.flux_type='intermolecular';
b.chem.flux_rate=[0 11;13 17]; b.chem.flux_type='intermolecular';
[merged_sys,merged]=merge_inp({part,part},{a,b});
result=test_close(result,'merged flux rate matrices',merged.chem.flux_rate,...
                  blkdiag(a.chem.flux_rate,b.chem.flux_rate),0,0,...
                  'original signed real rates and diagonal entries remain unchanged');
result=test_true(result,'merged common flux type',...
                 strcmp(merged.chem.flux_type,'intermolecular')&&numel(merged_sys.isotopes)==4,...
                 'matching flux types are accepted without converting to reactions');
for kind=1:3
    bad=b; field='flux_type';
    if kind==1, bad.chem.flux_type='intramolecular'; end
    if kind==2, bad.chem=rmfield(bad.chem,'flux_type'); end
    if kind==3
        bad.chem=rmfield(bad.chem,'flux_rate'); field='flux_rate';
    end
    rejected=false;
    try
        merge_inp({part,part},{a,bad});
    catch err
        rejected=contains(err.message,field);
    end
    result=test_true(result,['incompatible flux merge ' num2str(kind)],rejected,...
                     'mismatched or missing flux fields retain their specific native merge rejection');
end

% Cross-substance transfer invalidates independent-substance symmetry projectors
sys.isotopes={'1H','1H'}; inter=struct(); inter.zeeman.scalar={0,3};
inter.chem.parts={1,2}; inter.chem.concs=[1 1];
inter.chem.flux_rate=[0 5;7 0]; inter.chem.flux_type='intermolecular';
bas.approximation={'none','none'};
s=test_spin_system(sys,inter,bas);
F=relaxation(s); rho=state(s,'Lz',1); time=0.02;
current=evolution(s,1i*F,[],rho,time,1,'final');
result=test_close(result,'reduced cross-substance flux propagation',current,...
                  expm(full(F)*time)*rho,1e-12,0,...
                  'trajectory reduction retains source and destination in one coupled space');

% The unrelated cross-substance prohibition remains when no coupling is specified
bare=s; bare.chem.flux_rate=[]; bare.chem.flux_type=''; rejected=false;
try
    reduce(bare,1i*F,rho);
catch err
    rejected=strcmp(err.identifier,'Spinach:reduce:crossSubstanceGenerator');
end
result=test_true(result,'uncoupled substance guard retained',rejected,...
                 'only declared reactions or flux permit cross-substance generator blocks');

end



% Only actual transfer entries activate the flux consumer and reduction coupling
function result=inactive_flux_checks(result)

% No-op flux is allowed in either Zeeman formalism, as in the historical consumer
sys.magnet=14.1; sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,3};
inter.chem.flux_rate=zeros(2); inter.chem.flux_type='intermolecular';
bas.approximation={'none'};
for formalism={'zeeman-liouv','zeeman-hilb'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas);
    result=test_close(result,['zero flux in ' formalism{1}],...
                      relaxation(s),sparse(s.bas.offsets(end),s.bas.offsets(end)),0,0,...
                      'all-zero paired rates require no flux-specific formalism or basis');
end

% Caller-selected physical models must dispatch explicitly, not silently fall back
bas.formalism='sphten-liouv'; inter.chem.flux_rate=[0 3;7 0];
s=test_spin_system(sys,inter,bas);
for mode={'intermoleculr',''}
    s.chem.flux_type=mode{1}; rejected=false;
    try
        relaxation(s);
    catch err
        rejected=contains(err.message,'unrecognized type of magnetisation flux');
    end
    result=test_true(result,['invalid caller-selected flux model ' mode{1}],rejected,...
                     'the original physical dispatch rejects unsupported active-flux models');
end

% Zero rate matrices retain the same symmetry projectors as absent flux
inter.chem.parts={1,2}; inter.chem.concs=[1 1]; inter.chem.flux_rate=zeros(2);
bas.approximation={'none','none'}; sys.disable={'pt'};
s=test_spin_system(sys,inter,bas); bare=s;
bare.chem.flux_rate=[]; bare.chem.flux_type='';
rho=state(s,'Lz','1H'); dim=s.bas.offsets(end);
L=1i*spdiags((1:dim)',0,dim,dim);
result=test_true(result,'inactive flux preserves symmetry projectors',...
                 isequal(reduce(s,L,rho),reduce(bare,L,rho)),...
                 'an all-zero matrix does not couple independent substances');
L(7,3)=1; rejected=false;
try
    reduce(s,L,rho);
catch err
    rejected=strcmp(err.identifier,'Spinach:reduce:crossSubstanceGenerator');
end
result=test_true(result,'inactive flux preserves cross-substance guard',rejected,...
                 'all-zero rates do not authorize an undeclared cross-substance generator');

% Reject non-finite rates before they can enter a physical generator
sys.output='hush'; sys.parallel={'local',1};
invalid=inter; invalid.chem.flux_rate=ones(2,2,2); rejected=false;
try
    create(sys,invalid);
catch err
    rejected=contains(err.message,'square matrix');
end
result=test_true(result,'multidimensional flux rejected',rejected,...
                 'point-to-point rates must form a two-dimensional matrix');

% Reject non-finite rates
for rate=[NaN Inf -Inf]
    invalid=inter; invalid.chem.flux_rate(1,2)=rate; rejected=false;
    try
        create(sys,invalid);
    catch err
        rejected=contains(err.message,'finite');
    end
    result=test_true(result,['non-finite flux rate ' num2str(rate)],rejected,...
                     'every supplied point-to-point rate must be finite');
end

end


