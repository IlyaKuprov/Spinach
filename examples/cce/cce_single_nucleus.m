% Hahn ESEEM of one spin-half nucleus coupled to a spin-one centre.
% The conditional 0/+1 Hamiltonians give the analytic modulation in
% Eq. (29) of https://doi.org/10.1002/adts.202100254. Frequencies below
% are chosen model parameters, not a material-specific fit. Syntax:
%
%                     result=cce_single_nucleus()
%
% ilya.kuprov@weizmann.ac.il

function result=cce_single_nucleus()

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E3','13C'};
sys.output='hush'; sys.parallel={'local',1};

% Nuclear Zeeman frequency and conditional hyperfine components, rad/s
omega=-spin('13C')*sys.magnet; apar=0.3*abs(omega);
aper=0.6*abs(omega);

% Standard Spinach tensors: electron g, nuclear ppm, and coupling Hz
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=[0 0 aper;0 0 0;0 0 apar]/(2*pi);

% Compile metadata without constructing a full-bath basis
spin_system=create(sys,inter);

% Conditional coherence; static explicit operators include their angle
parameters.central_spins=1; parameters.bath_spins=2;
parameters.expansion='standard'; parameters.method='cce'; parameters.quantity='coherence'; parameters.order=1;
parameters.central.alpha=[0;1;0]; parameters.central.beta=[1;0;0];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.meanfield='none';
parameters.graph=struct('mode','complete');
parameters.times=linspace(0,8*pi/abs(omega),101);
parameters.sequence=struct('type','ideal','target',1,...
                           'operator',pi*[0 1 0;1 0 0;0 0 0]/2,...
                           'angle',1,'placement','fraction','at',0.5);
parameters.bath_state=struct('mode','mixed_product');
parameters.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^28);

% Run CCE and calculate an independent analytic echo
result=cce(spin_system,parameters);
omega1=sqrt((omega+apar)^2+aper^2);
reference=1-2*(aper/omega1)^2*sin(omega*parameters.times(:)/4).^2.*...
          sin(omega1*parameters.times(:)/4).^2;
result.reference=reference; result.reference_error=max(abs(result.value-reference));
fprintf('Single-nucleus analytic Hahn maximum error: %.3g\n',result.reference_error);

% Plot the complex CCE trace and the analytic real envelope
kfigure(); plot(1e6*result.times,real(result.value),'b-',...
                1e6*result.times,reference,'r--');
xlabel('total echo time / microseconds'); ylabel('normalised coherence');
legend('CCE1','analytic Eq. (29)'); kgrid();

end


