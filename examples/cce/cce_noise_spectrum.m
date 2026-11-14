% Ordered noise correlations and a finite-window spectrum of one nucleus.
% A precessing transverse spin has C(t)=cos(omega*t)/4 in a mixed state.
% This example compares that identity, not a phenomenological T1 fit.
% Syntax: result=cce_noise_spectrum()
%
% ilya.kuprov@weizmann.ac.il

function result=cce_noise_spectrum()

% Magnet and isotopes
sys.magnet=0.01; sys.isotopes={'E','13C'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);
spin_system=create(sys,inter);

% Additive connected noise of the transverse nuclear operator
x=[0 1;1 0]/2; omega=-spin('13C')*sys.magnet;
parameters.central_spins=1; parameters.bath_spins=2;
parameters.method='noise'; parameters.expansion='standard';
parameters.quantity='autocorrelation'; parameters.order=1;
parameters.central.rho=eye(2)/2; parameters.assumptions='labframe';
parameters.orientation=[0 0 0]; parameters.interaction_cutoff=Inf;
parameters.graph=struct('mode','complete'); parameters.meanfield='none';
parameters.times=linspace(0,10*2*pi/abs(omega),401);
parameters.sequence=struct([]); parameters.bath_state=struct('mode','mixed_product');
parameters.noise=struct('ordering','ordered','connected',true,...
                       'terms',struct('support',2,'operator',x));
parameters.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^28);

% Compare the additive expansion with an independent analytic correlation
result=cce(spin_system,parameters); reference=cos(omega*result.times(:))/4;
result.reference_error=max(abs(result.value-reference));
fprintf('Noise analytic cosine maximum error: %.3g.\n',result.reference_error);

% Transform the stated finite observation window, without a Markov assumption
parameters.spectrum=struct('frequency_grid',linspace(-2*abs(omega),2*abs(omega),401),...
    'convention','exp(+i*omega*t)','sidedness','two-sided','window','hann',...
    'time_domain','positive','weights','trapezoid');
spectrum=cce_spectrum(result.times,result.value,parameters); result.spectrum=spectrum;

% Plot the physical correlation and its finite-window transform
kfigure(); subplot(1,2,1); plot(1e6*result.times,real(result.value));
xlabel('lag / microseconds'); ylabel('C(t)'); kgrid();
subplot(1,2,2); plot(spectrum.omega/(2*pi*1e3),real(spectrum.value));
xlabel('frequency / kHz'); ylabel('finite-window spectrum / s'); kgrid();

end

