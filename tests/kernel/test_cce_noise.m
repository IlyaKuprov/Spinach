% Additive noise expansion against dense quantum and classical references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_noise()

% Prepare two bath spins and a nontracial product state
result=new_test_result('kernel/cce_noise','Additive bath noise',...
                       'Full order, disconnected cumulants, ordering and commuting limit.');
sys.isotopes={'1H','1H','1H'}; sys.magnet=0; sys.output='hush';
inter.zeeman.scalar={0,0,0}; inter.coupling.scalar=cell(3); s=create(sys,inter);
s.inter.zeeman.matrix{2}(3,3)=1.3; s.inter.zeeman.matrix{3}(3,3)=0.7;
s.inter.coupling.matrix{2,3}=diag([0.4 0.5 0.9]);
p=struct('central_spins',1,'bath_spins',[2 3],'method','noise',...
         'quantity','autocorrelation','expansion','standard','order',2,...
         'times',linspace(0,2,31),'assumptions','labframe',...
         'orientation',[0 0 0],'interaction_cutoff',Inf,'meanfield','none');
p.graph=struct('mode','complete'); p.bath_state=struct('mode','mixed_product');
p.sequence=struct([]); p.central=struct('alpha',[1;0],'beta',[0;1]);
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',1e8);
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
p.noise=struct('ordering','ordered','connected',true,...
               'terms',struct('support',{2,3},'operator',{x,2*z}));
r=cce(s,p); B=kron(x,eye(2))+2*kron(eye(2),z);
H=1.3*kron(z,eye(2))+0.7*kron(eye(2),z)+...
  0.4*kron(x,x)+0.5*kron(y,y)+0.9*kron(z,z);
reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    U=expm(-1i*H*p.times(k)); reference(k)=trace(U'*B*U*B)/4;
end
result=test_close(result,'C7 numerical check',r.value,reference,1e-12,1e-12,'Full-order additive exact reference.');
exact_error=max(abs(r.value-reference));

% Check factorization into independent connected single-spin correlations
s.inter.coupling.matrix{2,3}=zeros(3); p.order=1; first=cce(s,p);
p.order=2; second=cce(s,p);
reference=cos(1.3*p.times(:))/4+1;
result=test_close(result,'C7 numerical check',first.value,reference,1e-12,1e-12,'Disconnected independent sum.');
result=test_close(result,'C7 numerical check',second.value,first.value,1e-12,1e-12,'Disconnected pair cumulant vanishes.');

% Exercise unequal-time quantum ordering in a polarized bath
p.bath_state=struct('mode','explicit_product','rho',{{diag([0.8 0.2]),eye(2)/2}});
p.noise.terms=struct('support',2,'operator',x); ordered=cce(s,p);
reference=(cos(1.3*p.times(:))+0.6i*sin(1.3*p.times(:)))/4;
result=test_close(result,'C7 numerical check',ordered.value,reference,1e-12,1e-12,'Ordered complex quantum correlation.');
p.noise.ordering='symmetrized'; symmetric=cce(s,p);
result=test_close(result,'C7 numerical check',symmetric.value,real(reference),1e-12,1e-12,'Explicit anticommutator ordering.');

% Recover classical static variance for a commuting diagonal observable
p.noise.terms=struct('support',2,'operator',z); classical=cce(s,p);
result=test_close(result,'C7 numerical check',classical.value,repmat(0.16,numel(p.times),1),1e-12,1e-12,...
                  'Commuting classical connected variance.');
fprintf('Noise full-order error %.3g; quantum ordering error %.3g\n',...
        exact_error,max(abs(ordered.value-reference)));

end


