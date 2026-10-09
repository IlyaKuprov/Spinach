% Tests dense Zeeman chemistry and Hilbert matrix actions (T17 and T18).
% Syntax:
%
%                     result=test_cwdm_zeeman_kin()
%
% Outputs:
%
%    result - algebraic transport, trace, IME, and named rejection checks
%
% References are physical Hilbert matrices, tensor products, partial traces,
% and explicit singlet projectors. The tests do not obtain their expected
% reaction maps from spherical-tensor descriptors. Matrix comparisons use
% a 1e-12 absolute whole-array bound; integrated exchange uses 1e-10.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_zeeman_kin()

% Announce independent matrix-reference chemistry checks
result=new_test_result('kernel/cwdm_zeeman_kin','CWDM Zeeman chemistry T17-T18',...
                       'Dense reaction maps agree with physical matrix transport.');

% Build a permuted two-reactant product and a spin-free pool
sys.magnet=14.1; sys.isotopes={'1H','1H','1H','1H'};
inter.chem.parts={1,2,3:4,[]}; inter.chem.concs=[0.4 0.6 0 0];
inter.chem.reactions={struct('reactants',[1 2],'products',3,...
                             'matching',[1 4;2 3],'rate',7,'closure','product')};
bas.formalism='zeeman-liouv'; bas.approximation={'none','none','none','none'};
z=test_spin_system(sys,inter,bas);
rho_a=0.4*[0.7 0.1i;-0.1i 0.3]; rho_b=0.6*[0.2 0.05;0.05 0.8];
eta=hilb2liouv({rho_a,rho_b,zeros(4),0},'statevec');
K=kinetics(z); deriv=K(0,eta)*eta;
expected=hilb2liouv({-7*trace(rho_b)*rho_a,-7*trace(rho_a)*rho_b,...
                    7*kron(rho_b,rho_a),0},'statevec');
result=test_close(result,'product closure with permutation',deriv,expected,1e-12,0,...
                  'matching reverses the reactant spin order and preserves complex coherences');
concs=chem_concs(z,deriv); balance=[concs(1)+concs(3);concs(2)+concs(3)];
result=test_close(result,'product atom equivalents',balance,zeros(2,1),1e-12,0,...
                  'each reactant atom is transported into the product');
fprintf('CWDM_ZEEMAN_PRODUCT error=%.16g atom_balance=%.16g\n',norm(deriv-expected),norm(balance));

% Evaluate concentrations from the propagated state rather than initial inputs
result=test_close(result,'Zeeman trace extraction',chem_concs(z,eta),[0.4 0.6 0 0],1e-14,0,...
                  'trace functionals replace spherical-tensor unit coordinates');
spatial=[eta;2*eta]; spatial_k=K(0,spatial);
result=test_close(result,'spatial product closure',spatial_k*spatial,[expected;4*expected],1e-12,0,...
                  'each voxel gets its own quadratic mass-action extent');
zero_eta=hilb2liouv({zeros(2),rho_b,zeros(4),0},'statevec');
result=test_close(result,'zero-population product',K(0,zero_eta)*zero_eta,0*eta,1e-12,0,...
                  'an empty reactant gives zero extent without concentration division');

% Integrate the product reaction through the production state-dependent stepper
options=odeset('RelTol',1e-12,'AbsTol',1e-14);
[~,reference]=ode45(@product_rhs,[0 0.01],eta,options);
final=eta; dt=0.001;
for n=1:10
    final=step(z,{@(t,x)1i*K(t,x),(n-1)*dt,'RKMK4'},final,dt);
end
result=test_close(result,'Zeeman nonlinear integration',final,reference(end,:).',1e-10,0,...
                  'RKMK4 agrees with the independent Hilbert product-state ODE');
fprintf('CWDM_ZEEMAN_RKMK error=%.16g\n',norm(final-reference(end,:).'));

% Compare additive arrival with the tensor-product expansion without cross orders
z.chem.reactions{1}.closure='additive'; K=kinetics(z);
arrival=7*(trace(rho_b)*kron(eye(2)/2,rho_a)+...
           trace(rho_a)*kron(rho_b,eye(2)/2)-...
           trace(rho_a)*trace(rho_b)*eye(4)/4);
expected=hilb2liouv({-7*trace(rho_b)*rho_a,-7*trace(rho_a)*rho_b,arrival,0},'statevec');
result=test_close(result,'additive matrix reference',K(0,eta)*eta,expected,1e-12,0,...
                  'the identity arrives once and cross-reactant polarisation products are omitted');
z.chem.reactions{1}.rate=@(t)7*(1+t); K=kinetics(z);
result=test_close(result,'time-dependent dense map',K(0.3,eta)*eta,1.3*expected,1e-12,0,...
                  'time rates use the same shared chemical generator');

% Reject every Hilbert mass-action closure by its complete documented message
bas.formalism='zeeman-hilb'; h=test_spin_system(sys,inter,bas);
for closure={'additive','product'}
    h.chem.reactions{1}.closure=closure{1}; rejected=false;
    try
        kinetics(h);
    catch err
        rejected=strcmp(err.identifier,'Spinach:kinetics:hilbertMassAction')&&...
                 strcmp(err.message,'mass-action chemistry is not supported in zeeman-hilb formalism.');
    end
    result=test_true(result,['T18 Hilbert ' closure{1}],rejected,...
                     'matrix-valued mass-action propagation is rejected rather than approximated');
end

% Test partial trace and unpolarised spin arrival with unequal dimensions
inter.chem.reactions={struct('reactants',3,'products',1,'matching',[3 1],'rate',2),...
                      struct('reactants',1,'products',3,'matching',[1 3],'rate',3)};
singlet=[0;1;-1;0]/sqrt(2); rho_c=0.2*(singlet*singlet');
expected_a=2*0.2*eye(2)/2-3*rho_a;
expected_c=-2*rho_c+3*kron(rho_a,eye(2)/2);
for formalism={'zeeman-liouv','zeeman-hilb'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas); K=kinetics(s);
    if strcmp(formalism{1},'zeeman-liouv')
        rho=hilb2liouv({rho_a,rho_b,rho_c,0},'statevec');
        deriv=K*rho; expected=hilb2liouv({expected_a,zeros(2),expected_c,0},'statevec');
    else
        rho=blkdiag(rho_a,rho_b,rho_c,0);
        deriv=K(0,rho); expected=blkdiag(expected_a,zeros(2),expected_c,0);
    end
    result=test_close(result,['unequal transport ' formalism{1}],deriv,expected,1e-12,0,...
                      'unmatched source spins are traced out and unmatched product spins arrive unpolarised');
end

% Evolve first-order Hilbert exchange with an independent matrix-exponential reference
sys.isotopes={'1H','1H'}; inter=struct('temperature',298);
inter.zeeman.scalar={1,2}; inter.chem.parts={1,2}; inter.chem.concs=[0.4 0.6];
inter.chem.reactions={struct('reactants',1,'products',2,'matching',[1 2],'rate',20),...
                      struct('reactants',2,'products',1,'matching',[2 1],'rate',10)};
bas.formalism='zeeman-hilb'; bas.approximation={'none','none'};
h=test_spin_system(sys,inter,bas); K=kinetics(h); rho=blkdiag(rho_a,rho_b);
deriv=K(0,rho); trace_error=abs(trace(deriv));
result=test_close(result,'T17 production exchange',deriv,...
                  blkdiag(-20*rho_a+10*rho_b,20*rho_a-10*rho_b),1e-12,0,...
                  'Hilbert kinetics returns a matrix derivative, not a conjugation Hamiltonian');
result=test_close(result,'T17 production traces',chem_concs(h,rho),[0.4 0.6],1e-14,0,...
                  'Hilbert concentrations are the local traces');
options=odeset('RelTol',1e-12,'AbsTol',1e-14);
[~,trajectory]=ode45(@(t,x)reshape(K(t,reshape(x,4,4)),[],1),[0 0.07],rho(:),options);
final=reshape(trajectory(end,:),4,4);
expected=expm(kron([-20 10;20 -10],eye(4))*0.07)*hilb2liouv({rho_a,rho_b},'statevec');
result=test_close(result,'integrated Hilbert exchange',...
                  hilb2liouv({final(1:2,1:2),final(3:4,3:4)},'statevec'),expected,1e-10,0,...
                  'the matrix RHS integrates to the full-component exchange exponential');
fprintf('CWDM_T17_PRODUCTION trace=%.16g\n',trace_error);

% Retain time-dependent first-order rates in the Hilbert matrix interface
h.chem.reactions{1}.rate=@(t)20*(1+t); K=kinetics(h);
result=test_close(result,'Hilbert time-dependent exchange',K(0.2,rho),...
                  blkdiag(-24*rho_a+10*rho_b,24*rho_a-10*rho_b),1e-12,0,...
                  'the matrix callback evaluates the rate at the requested time');

% Apply Hilbert IME through local trace-one targets and trace-preserving maps
unit_h=h; unit_h.chem.concs(:)=1; target=equilibrium(unit_h);
unit=reshape(eye(2),4,1); bare=blkdiag(-3*(eye(4)-unit*unit'/2),-4*(eye(4)-unit*unit'/2));
R=thermalize(h,bare,[],[],target,'IME'); stationary=R(0,equilibrium(h));
result=test_close(result,'Hilbert IME stationary',stationary,zeros(4),1e-12,0,...
                  'the matrix target follows each block concentration');
deriv=R(0,rho);
expected=blkdiag(-3*(rho_a-trace(rho_a)*target(1:2,1:2)),...
                 -4*(rho_b-trace(rho_b)*target(3:4,3:4)));
result=test_close(result,'Hilbert IME action',deriv,expected,1e-12,0,...
                  'matrix IME implements R times the deviation from concentration-weighted equilibrium');
result=test_close(result,'Hilbert IME trace',chem_concs(h,deriv),[0 0],1e-12,0,...
                  'relaxation does not change concentrations');
fprintf('CWDM_HILBERT_IME stationary=%.16g trace=%.16g\n',norm(stationary,'fro'),norm(chem_concs(h,deriv)));

% Reject matrix coherences between species rather than silently projecting them out
bad=rho; bad(1,3)=1;
calls={@()K(0,bad),@()chem_concs(h,bad)};
ids={'Spinach:hilb_action:crossSubstance','Spinach:chem_concs:crossSubstance'};
for n=1:numel(calls)
    rejected=false;
    try
        calls{n}();
    catch err
        rejected=strcmp(err.identifier,ids{n});
    end
    result=test_true(result,ids{n},rejected,'unphysical off-diagonal species blocks are not discarded');
end
rejected=false;
try
    thermalize(h,bare,[],[],equilibrium(h),'IME');
catch err
    rejected=strcmp(err.identifier,'Spinach:thermalize:targetConcentration');
end
result=test_true(result,'Hilbert weighted IME target',rejected,...
                 'matrix IME requires unit-concentration targets, not doubly weighted targets');

% Verify named and user-supplied spin-selective loss and spin-free arrival
sys.isotopes={'E','E'}; inter=struct(); inter.chem.parts={1:2,[]}; inter.chem.concs=[1 0];
inter.chem.reactions={struct('reactants',1,'products',2,'matching',zeros(0,2),...
                             'rate',7,'selector',{{'singlet',[1 2]}})};
P=singlet*singlet'; density=eye(4)/4;
expected=blkdiag(-7*(P*density+density*P)/2,7*trace(P*density*P));
for formalism={'zeeman-liouv','zeeman-hilb'}
    bas.formalism=formalism{1}; s=test_spin_system(sys,inter,bas);
    for user_selector=[false true]
        if user_selector
            s.chem.reactions{1}.selector={hilb2liouv(P,'left'),hilb2liouv(P,'right')};
        end
        K=kinetics(s);
        if strcmp(formalism{1},'zeeman-liouv')
            deriv=K*hilb2liouv({density,0},'statevec');
            reference=hilb2liouv({expected(1:4,1:4),expected(5,5)},'statevec');
        else
            deriv=K(0,blkdiag(density,0)); reference=expected;
        end
        result=test_close(result,[formalism{1} ' selector ' num2str(user_selector)],deriv,reference,1e-12,0,...
                          'Haberkorn loss feeds the traced singlet population to the spin-free product');
    end
end
fprintf('CWDM_T18_COMPLETE failures=%d\n',numel(result.failures));
end

% Independent physical matrix RHS for the permuted bimolecular product
function deriv=product_rhs(~,eta)

% Read the two reacting physical density matrices
rho_a=reshape(eta(1:4),2,2); rho_b=reshape(eta(5:8),2,2);

% Form mass-action losses and the reversed-order tensor-product arrival
dot_a=-7*trace(rho_b)*rho_a; dot_b=-7*trace(rho_a)*rho_b;
dot_c=7*kron(rho_b,rho_a);
deriv=[dot_a(:);dot_b(:);dot_c(:);0];

end


