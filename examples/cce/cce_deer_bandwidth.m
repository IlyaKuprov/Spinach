% Hahn echo with a frequency-selected bath pump and finite bandwidth.
% Uses cce_controls explicit Hz site offsets; no inferred resonance model.
% The deliberately small model is not a published material reproduction.
% Syntax: result=cce_deer_bandwidth()
%
% ilya.kuprov@weizmann.ac.il

function result=cce_deer_bandwidth()

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','E','E'};
sys.output='hush'; sys.parallel={'local',1};

% Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
inter.coordinates={[0 0 0];[25 0 0];[0 35 0]};
inter.zeeman.matrix=repmat({zeros(3)},1,3);
inter.coupling.matrix=cell(3);

% Explicit illustrative hyperfine and bath coupling tensors in Hz
inter.coupling.matrix{1,2}=diag([0 0 20000]);
inter.coupling.matrix{1,3}=diag([0 0 12000]);
inter.coupling.matrix{2,3}=diag([200 200 -400]);

% Coordinates define the toy geometry; all couplings above are explicit
inter=rmfield(inter,'coordinates');

% Compile physical metadata without a full-bath basis
spin_system=create(sys,inter);

% Conditional central-spin Hahn coherence and a complete small-bath graph
parameters.central_spins=1; parameters.bath_spins=2:3;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='coherence'; parameters.order=2;
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,200e-6,41);
parameters.sequence=struct('type','ideal','target',1,'axis','x',...
                           'angle',pi,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Select the first bath resonance with a rectangular pump window
parameters.sequence(2)=struct('type','ideal','target',2,'axis','x',...
                              'angle',pi,'placement','fraction','at',0.25);
parameters.controls.frequencies=struct('sites',[2 3],'values',[0 1e6],'units','Hz');
parameters.controls.pulses=struct('event',2,'scope',[2 3],'centre',0,...
    'bandwidth',2e5,'offset',0,'units','Hz','shape','rectangular');

% Run the cluster expansion and report its physical input boundary
result=cce(spin_system,parameters);
fprintf('cce_deer_bandwidth: %d/%d valid values; finite-cluster demonstration.\n',...
        nnz(result.valid),numel(result.valid));
fprintf('Pump centre 0 Hz, bandwidth 200 kHz; offsets [0,1 MHz].\n');

% Plot the conditional complex coherence
kfigure(); plot(1e6*result.times,real(result.value),'b-',...
                1e6*result.times,imag(result.value),'r--');
xlabel('total echo time / microseconds'); ylabel('normalised coherence');
legend('real','imaginary'); kgrid();

end

