% Selector loss and CIDNP product transport tests T14-T15. Syntax:
%
%                   result=test_cwdm_selectors()
%
% Outputs:
%
%    result - comparisons with independent Hilbert-space projector algebra
%
% The radical pair contains two electrons and a nucleus. Singlet-channel
% products retain the nucleus; triplet products are spin-free. Physical
% Hilbert matrices use the converter's trace/D convention explicitly.
% T15 includes noncommuting electron-nuclear evolution and an augmented
% matrix exponential for the integrated projected product source.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_selectors()

% Build independently labelled radical-pair and diamagnetic product blocks
fprintf('TESTING: CWDM selectors and CIDNP transport T14-T15\n');
result=new_test_result('kernel/cwdm_selectors','CWDM selectors T14-T15',...
                      'Selective reaction loss and nuclear product arrival against Hilbert algebra.');
sys.magnet=1; sys.isotopes={'E','E','1H','1H'};
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.chem.parts={1:3,4,[]}; inter.chem.concs=[0.7 0 0];
singlet_rx=struct('reactants',1,'products',2,'matching',[3 4],'rate',2,...
                  'selector',{{'singlet',[1 2]}});
triplet_rx=struct('reactants',1,'products',3,'matching',zeros(0,2),'rate',3,...
                  'selector',{{'triplet',[1 2]}});
inter.chem.reactions={singlet_rx,triplet_rx};
bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};
s=basis(create(sys,inter),bas); K=kinetics(s); P=full(sphten2zeeman(s));
source=1:64; product=65:68; sink=69;

% Form physical projectors without Spinach operator construction
singlet_ket=[0;1;-1;0]/sqrt(2);
singlet_proj=kron(singlet_ket*singlet_ket',eye(2)); triplet_proj=eye(8)-singlet_proj;
rng(14); trial=randn(8)+1i*randn(8);
rho=trial*trial'; rho=0.7*rho/trace(rho);
eta=zeros(69,1); eta(source)=P(source,source)\(8*rho(:));

% Compare selective loss and arrival with direct Hilbert multiplication
selected=singlet_proj*rho*singlet_proj; nuclear=zeros(2);
for n=1:4
    block=(2*n-1):(2*n); nuclear=nuclear+selected(block,block);
end
expected=-singlet_rx.rate*(singlet_proj*rho+rho*singlet_proj)/2 ...
         -triplet_rx.rate*(triplet_proj*rho+rho*triplet_proj)/2;
deriv=K*eta; observed=reshape(P(source,source)*deriv(source),8,8)/8;
result=test_close(result,'T14 Haberkorn drain',observed,expected,1e-12,0,...
                  'selective loss is the anticommutator with each Hilbert projector');
result=test_close(result,'T15 projected partial trace',...
                  reshape(P(product,product)*deriv(product),2,2)/2,...
                  singlet_rx.rate*nuclear,1e-12,0,...
                  'electron tracing retains the projected nuclear density matrix');
expected=[-2*trace(singlet_proj*rho)-3*trace(triplet_proj*rho),2*trace(singlet_proj*rho),3*trace(triplet_proj*rho)];
result=test_close(result,'selective concentration balance',chem_concs(s,deriv),expected,1e-12,0,...
                  'channel populations, rather than unselected mass-action concentrations, determine rates');
result=test_close(result,'tracked selective conservation',sum(chem_concs(s,deriv)),0,1e-12,0,...
                  'every lost radical pair arrives in one tracked product');

% Reject stationary tracked products in the full-space RYDMR resolvent
parameters.tol=1e-12; zero=sparse(size(K,1),size(K,2)); rejected=false;
try
    rydmr(s,parameters,zero,zero,K);
catch err
    rejected=strcmp(err.identifier,'Spinach:rydmr:trackedProducts');
end
result=test_true(result,'RYDMR tracked product guard',rejected,...
                 'tracked populations require time-domain propagation rather than a stationary full-space source');

% Retain the analytic singlet yield for an untracked recombination sink
loss_system=s;
for n=1:numel(loss_system.chem.reactions)
    loss_system.chem.reactions{n}.products=[];
    loss_system.chem.reactions{n}.matching=zeros(0,2);
end
yield=rydmr(loss_system,parameters,zero,zero,kinetics(loss_system));
result=test_close(result,'RYDMR untracked singlet yield',yield,1,1e-12,0,...
                  'without spin mixing the prepared singlet recombines entirely through its channel');

% Compare Jones-Hore loss on a density with singlet-triplet coherences
jones=s; jones.chem.reactions{1}.selector{1}='jones-hore-singlet';
jones.chem.reactions{2}.selector{1}='jones-hore-triplet';
jones_gen=kinetics(jones); deriv=jones_gen*eta;
expected=-2*(rho-triplet_proj*rho*triplet_proj)-3*(rho-singlet_proj*rho*singlet_proj);
result=test_close(result,'T14 Jones-Hore drain',...
                  reshape(P(source,source)*deriv(source),8,8)/8,expected,1e-12,0,...
                  'each channel retains the complementary projected state');
result=test_close(result,'Jones-Hore projected arrival',deriv(product),...
                  K(product,:)*eta,1e-12,0,'channel dephasing changes loss, not projected arrival');

% The user-specified local projector pair reproduces the named channel
custom=s; left=kron(eye(8),singlet_proj); right=kron(singlet_proj.',eye(8));
custom.chem.reactions{1}.selector={P(source,source)\(left*P(source,source)),...
                                  P(source,source)\(right*P(source,source))};
result=test_close(result,'user projector pair',kinetics(custom),K,1e-12,0,...
                  'local left/right product matrices implement the named singlet selector');

% Integrate a noncommuting hyperfine Hamiltonian and projected source exactly
spin=pauli(2);
hilb_ham=0.2*kron(kron(spin.z,eye(2)),spin.z)+0.3*kron(kron(spin.x,eye(2)),spin.x);
H=0.2*operator(s,{'Lz','Lz'},{1,3})+0.3*operator(s,{'Lx','Lx'},{1,3});
loss=-1i*(kron(eye(8),hilb_ham)-kron(hilb_ham.',eye(8)))...
     -(2*(kron(eye(8),singlet_proj)+kron(singlet_proj.',eye(8)))...
       +3*(kron(eye(8),triplet_proj)+kron(triplet_proj.',eye(8))))/2;
arrival=zeros(5,64);
for n=1:64
    element=zeros(8); element(n)=1; selected=singlet_proj*element*singlet_proj;
    nuclear=zeros(2);
    for k=1:4
        block=(2*k-1):(2*k); nuclear=nuclear+selected(block,block);
    end
    arrival(:,n)=[2*nuclear(:);3*trace(triplet_proj*element)];
end
reference=expm([loss zeros(64,5);arrival zeros(5)]*0.43)*[rho(:);zeros(5,1)];
final=expm(full(K-1i*H)*0.43)*eta;
observed=[P(source,source)*final(source)/8;P(product,product)*final(product)/2;final(sink)];
result=test_close(result,'T15 integrated product source',observed,reference,1e-10,0,...
                  'the augmented Hilbert exponential integrates kS times the projected nuclear source');

% Automatic reduction must preserve chemical links into empty product blocks
reduced=evolution(s,H+1i*K,[],eta,0.43,1,'final');
result=test_close(result,'T15 reduced propagation',reduced,final,1e-10,0,...
                  'substance-local spin irreps are not independent under chemical transport');
trimmed=s; trimmed.sys.enable=[trimmed.sys.enable {'zte'}];
reduced=evolution(trimmed,H+1i*K,[],eta,0.43,1,'final');
result=test_close(result,'T15 ZTE product arrival',reduced,final,1e-10,0,...
                  'full-generator reductions retain initially empty reachable product coordinates');

% Detect the independently integrated nuclear product polarisation
coil=coil_state(s,'Lz',4,'exact');
expected=trace(spin.z*reshape(reference(product),2,2));
result=test_close(result,'T15 nuclear polarisation',coil'*final,expected,1e-10,0,...
                  'the product coil detects the transported nuclear expectation value');
fprintf('CWDM_T15_HILBERT_ERROR %.12g POLARISATION_ERROR %.12g\n',...
        norm(observed-reference),abs(coil'*final-expected));

end


