% A central spin coupled to a phenomenologically relaxing spin-half bath.
% ME-gCCE is compared with an independent vectorised dense Liouvillian.
% The supplied rate is not a prediction of microscopic T1. See
% https://doi.org/10.1103/PhysRevLett.132.250401. Syntax:
%
%                       result=mecce_dissipative_bath()
%
% ilya.kuprov@weizmann.ac.il

function result=mecce_dissipative_bath()

% Standard Spinach system and anisotropic hyperfine tensor in Hz
sys.magnet=0; sys.isotopes={'E','13C'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([1100 900 2200]);
spin_system=create(sys,inter);

% Full density ME-gCCE with an explicitly supplied nuclear jump rate
parameters.central_spins=1; parameters.bath_spins=2;
parameters.method='megcce'; parameters.expansion='standard';
parameters.quantity='density'; parameters.order=1;
v=[1;1]/sqrt(2); parameters.central.rho=v*v';
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete'); parameters.times=linspace(0,200e-6,31);
parameters.sequence=struct([]); parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);
jump=[0 1;0 0]; rate=3000;
parameters.jumps=struct('support',2,'operator',jump,'rate',rate);

% Propagate the reduced density with the cluster implementation
result=cce(spin_system,parameters);

% Independently construct the column-vectorised dense master equation
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
h=2*pi*(1100*kron(x,x)+900*kron(y,y)+2200*kron(z,z));
a=kron(eye(2),jump); q=a'*a; unit=eye(4);
liouvillian=-1i*(kron(unit,h)-kron(h.',unit))+...
    rate*(kron(conj(a),a)-(kron(unit,q)+kron(q.',unit))/2);
rho0=kron(parameters.central.rho,eye(2)/2); reference=zeros(size(result.value));
for n=1:numel(result.times)

    % Dense exponential followed by an explicit nuclear partial trace
    rho=reshape(expm(liouvillian*result.times(n))*rho0(:),4,4); reduced=zeros(2);
    for k=1:2

        % Tensor order is central electron, then bath nucleus
        ids=k:2:4; reduced=reduced+rho(ids,ids);
    end
    reference(n,:)=reduced(:).';
end
result.reference_error=max(abs(result.value-reference),[],'all');
fprintf('ME-gCCE dense Liouvillian maximum error: %.3g; supplied rate %.1f /s.\n',...
        result.reference_error,rate);

% Display central coherence without fitting a microscopic relaxation time
kfigure(); plot(1e6*result.times,real(result.value(:,3)),'b-',...
                1e6*result.times,real(reference(:,3)),'r--');
xlabel('time / microseconds'); ylabel('central coherence');
legend('ME-gCCE1','dense Liouvillian'); kgrid();

end

