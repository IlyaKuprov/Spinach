% Spatial reaction transport with a spin-free solvent pool. Syntax:
%
%                    result=test_cwdm_spatial()
%
% Outputs:
%
%    result - checks of mass action, frozen-history propagation, and
%             additive spin transport on a two-cell model
%
% This compact test is not full-chip reacting_flow_nmr acceptance. The
% concentration reference uses the original asymmetric allocation of
% the product unit source; the kernel shares that source equally.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_spatial()

% Construct two bimolecular channels and an inert spin-free solvent
fprintf('TESTING: CWDM spatial chemistry and concentration history\n');
result=new_test_result('kernel/cwdm_spatial','CWDM spatial chemistry',...
                      'Two-cell chemistry and transport against independent mass action.');
sys.magnet=1; sys.isotopes={'1H','1H','1H','1H','1H','1H'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.chem.parts={1,2,3:4,5:6,[]}; inter.chem.concs=[0.6 0.5 0 0 18.1];
inter.chem.reactions={struct('reactants',[1 2],'products',3,...
                            'matching',[1 3;2 4],'rate',2),...
                      struct('reactants',[1 2],'products',4,...
                            'matching',[1 5;2 6],'rate',1)};
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none','none','none'};
s=basis(create(sys,inter),bas); K=kinetics(s);
chem=kill_spin(s,1:s.comp.nspins); K_chem=kinetics(chem);
result=test_true(result,'five scalar pools',isequal(chem.bas.nstates,ones(5,1)),...
                 'tracing every spin preserves all five chemical substances');

% Compare spatial unit derivatives with an independent mass-action formula
concs=[0.6 0.2;0.5 0.3;0 0.1;0 0.2;18.1 17];
GF=[-0.2 0.3;0.2 -0.3]; expected=zeros(5,2);
for n=1:2
    expected(:,n)=concs(1,n)*concs(2,n)*[-3;-3;2;1;0];
end
expected=expected+concs*GF.';
G=K_chem(0,concs(:))+kron(GF,speye(5));
result=test_close(result,'spatial mass action',reshape(G*concs(:),5,2),...
                  expected,1e-12,0,'reaction is local and transport acts on every pool');

% Quantify the difference between two frozen-rate unit-source allocations
reference=cell(2,1);
for n=1:2
    a=concs(1,n); b=concs(2,n);
    reference{n}=[-3*b 0 0 0 0;0 -3*a 0 0 0;...
                  0 2*a 0 0 0;0 a 0 0 0;0 0 0 0 0];
end
reference=blkdiag(reference{:})+kron(GF,speye(5));
result=test_close(result,'original instantaneous mass action',...
                  reference*concs(:),G*concs(:),1e-12,0,...
                  'both frozen generators agree on the state where they were assembled');
dt=0.02; old_next=expm(reference*dt)*concs(:);
new_next=step(chem,1i*G,concs(:),dt);
fprintf('CWDM_SPATIAL_FROZEN_DIFFERENCE %.12g\n',norm(new_next-old_next));
result=test_true(result,'distinct frozen trajectories',norm(new_next-old_next)>1e-8,...
                 'equal instantaneous derivatives do not imply equal frozen exponentials');
result=test_close(result,'kernel frozen expm',new_next,expm(G*dt)*concs(:),1e-12,0,...
                  'the production step applies the compiled frozen generator');

% Embed voxel concentrations and polarisation in the full direct sum
units=s.bas.offsets(1:end-1)+1; dim=s.bas.offsets(end);
eta=zeros(dim,2); eta(units,:)=concs;
eta=eta+0.1*coil_state(s,'Lz',1)*concs(1,:)...
       +0.2*coil_state(s,'Lz',2)*concs(2,:);
full_gen=K(0,eta(:))+kron(GF,speye(dim));
result=test_close(result,'full spin unit derivatives',chem_concs(s,full_gen*eta(:)),...
                  expected.',1e-12,0,'spin orders cannot alter nonselective mass action');
spin_next=step(s,1i*full_gen,eta(:),dt);
result=test_close(result,'full spin concentration history',chem_concs(s,spin_next),...
                  reshape(new_next,5,2).',1e-12,0,'full and traced systems share their unit dynamics');
result=test_close(result,'solvent no reaction',...
                  full_gen(units(5),:)*eta(:),expected(5,1),1e-12,0,...
                  'the solvent moves between cells without reacting');

end


