% Ordered noise correlations and a finite-window spectrum of one nucleus.
% A precessing transverse spin has C(t)=cos(omega*t)/4 in a mixed state.
% The regression test checks that identity, not a phenomenological T1 fit.
% Syntax:
%
%             [answer,diagnostics]=cce_noise_small()
%
% Returns the correlation column in answer, in squared spin-operator units,
% and the calculation diagnostics, including the spectrum, in diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_noise_small()

% Magnet and isotopes
sys.magnet=0.01;
sys.isotopes={'E','13C'};
sys.output='hush';
sys.parallel={'local',1};
sys.tols.cce_atol=1e-9;
sys.tols.cce_rtol=1e-9;
sys.tols.cce_mem=2^28;
inter.zeeman.matrix={zeros(3),zeros(3)};
inter.coupling.matrix=cell(2);

% Spinach housekeeping
spin_system=create(sys,inter);

% Additive connected noise of the transverse nuclear operator
omega=-spin('13C')*sys.magnet;
parameters.observed=1;
parameters.cce.order=1;
parameters.preparation.state='mixed';
parameters.times=linspace(0,10*2*pi/abs(omega),401);
parameters.noise.descriptors=struct('operators','Lx','spins',2,'coefficient',1);
parameters.noise.ordering='ordered';
parameters.noise.connected=true;
parameters.noise.manifold='bath_only';

% Extend the stationary mixed-state correlation without a Markov assumption
parameters.spectrum=struct('frequency_grid',...
                           linspace(-2*abs(omega),2*abs(omega),401),...
                           'convention','exp(+i*omega*t)',...
                           'sidedness','two-sided','window','hann',...
                           'time_domain','positive','weights','trapezoid',...
                           'stationary',true);
[answer,diagnostics]=clexpn(spin_system,@noise_correlation,parameters,'labframe');

end


