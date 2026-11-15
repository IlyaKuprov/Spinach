% Orientation-averaged molecular Hahn echo with explicit quadrature weights.
% Two deuterons illustrate tensor rotation; three axes are not converged powder.
% The deliberately small model is not a published material reproduction.
% Syntax: result=cce_molecular_powder()
%
% ilya.kuprov@weizmann.ac.il

function result=cce_molecular_powder()

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','2H','2H'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[1.45 0 0];[-0.725 1.2557 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=[0 0 0;0 0 0;40000 0 20000];
inter.coupling.matrix{1,3}=[0 0 0;0 0 0;20000 0 -10000];
inter.coupling.matrix{2,3}=diag([500 500 -1000]);
inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
inter.coupling.matrix{3,3}=diag([-1 -1 2])*5e4;

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Generalized observed-spin Hahn signal and a complete small-bath graph
parameters.observed_spins=1; parameters.bath_spins=2:3;
parameters.expansion='standard'; parameters.method='gcce'; parameters.quantity='signal'; parameters.order=2;
parameters.transition.alpha=[1;0]; parameters.transition.beta=[0;1];

% Rotated tensors mix observed states, so retain their full generalized dynamics
parameters.rho0=ones(2)/2; parameters.normalized=true;
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,41);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);


% Average explicit orientations with positive normalised quadrature weights
orientations=[0 0 0;0 pi/2 0;pi/2 pi/2 0]; weights=ones(3,1)/3;
values=complex(zeros(numel(parameters.times),3));
valid=true(size(values)); bounds=zeros(size(values)); raw={};
for n=1:3

    % Rotate all anisotropic tensors consistently before the CCE calculation
    parameters.orientation=orientations(n,:); current=cce(spin_system,parameters);
    values(:,n)=current.value; valid(:,n)=current.valid;
    bounds(:,n)=current.abs_error_bound; raw=[raw current.raw]; %#ok<AGROW>
    if n==1
        result=current; result.value(:)=0; result.abs_error_bound(:)=0;
        result.by_order(:)=0; result.by_order_error(:)=0;
        result.elapsed_seconds=0;
    end

    % Average response arrays and numerical bounds with the same weights
    result.value=result.value+weights(n)*current.value;
    result.abs_error_bound=result.abs_error_bound+weights(n)*current.abs_error_bound;
    result.by_order=result.by_order+weights(n)*current.by_order;
    result.by_order_error=result.by_order_error+weights(n)*current.by_order_error;
    result.valid=result.valid&current.valid;
    result.by_order_valid=result.by_order_valid&current.by_order_valid;
    failed=~current.by_order_valid;
    result.by_order_reason(failed)=current.by_order_reason(failed);
    result.elapsed_seconds=result.elapsed_seconds+current.elapsed_seconds;
end

% Keep orientation-resolved values and recompute the averaged diagnostics
result.response=result.value; result.reason=result.by_order_reason(:,:,end);
result.orientation_values=values;
result.orientation_abs_error_bound=bounds; result.orientation_weights=weights;
result.orientation_valid=valid; result.orientations=orientations;
result.raw=raw; result.budget_satisfied=result.valid&...
    result.abs_error_bound<=parameters.numerics.atol+parameters.numerics.rtol*abs(result.value);
model=cce_prepare(spin_system,parameters);
result.diagnostics=cce_diagnostics(model,result.raw,result,parameters);
fprintf('Three-axis powder illustration; not a converged material powder curve.\n');

% Plot the projected complex signal
kfigure(); plot(1e6*result.times,real(result.value),'b-',...
                1e6*result.times,imag(result.value),'r--');
xlabel('total echo time / microseconds'); ylabel('normalised echo amplitude E(2τ)');
legend('real','imaginary'); kgrid();

end


