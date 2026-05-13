% Tests voitlander() on an isotropic one-electron field-swept line. Syntax:
%
%                    result=test_dynamic_voitlander()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test uses an isotropic spin-half electron, for which all triangle
% vertices and subdivision midpoints have the same transition field.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_voitlander()

% Announce the test target
fprintf('TESTING: Voitlander spherical-triangle integration\n');

% State the Voitlander target of the test
result=new_test_result('kernel/dynamic_voitlander',...
                       'Adaptive Voitlander triangle integral',...
                       'voitlander() must integrate a constant isotropic transition as triangle area times a Lorentzian line.');

% Build an isotropic one-electron Hilbert-space system
sys.magnet=1;
sys.isotopes={'E'};
inter.zeeman.matrix={2.0023*eye(3)};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Set field-swept EPR parameters with a deliberately loose recursion gate
parameters.spins={'E'};
parameters.mw_freq=9.5e9;
parameters.fwhm=2e-3;
parameters.window=[0.33 0.35];
parameters.npoints=9;
parameters.tm_tol=0;
parameters.rspt_order=Inf;
parameters.int_tol=1e9;
parameters.pp_tol=(parameters.window(2)-parameters.window(1))/(2*(parameters.npoints-1));
parameters.orientation=[0 0 0];
parameters.rho0=-state(spin_system,'Lz','E');

% Get Zeeman, coupling, and microwave Hamiltonians
[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'));
[Iz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'));
Ic=(Ic+Ic')/2;
Iz=(Iz+Iz')/2;
Hmw=(state(spin_system,'L+','E')+state(spin_system,'L-','E'))/2;

% Find the isotropic transition at one orientation
[tf,tm,tw,pd]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw);
result=test_true(result,'voitlander transition count',isscalar(tf),...
                 'an isolated isotropic electron has one allowed EPR transition in the selected field window');

% Define the positive octant spherical triangle
r1=[1;0;0];
r2=[0;1;0];
r3=[0;0;1];

% Keep all field points inside the internal six-width support window
b_axis=tf+linspace(-2*tw,2*tw,parameters.npoints);

% Integrate the triangle using the production routine
spec=voitlander(spin_system,parameters,r1,r2,r3,...
                tf,tf,tf,tm,tm,tm,tw,tw,tw,pd,pd,pd,...
                Ic,Iz,Qc,Qz,Hmw,b_axis);

% Build the analytic Lorentzian reference for a constant transition
line_width=tw/2;
line_shape=(tm/(2*pi*line_width))./(1+((b_axis-tf)/line_width).^2);
reference=sphtarea(r1,r2,r3)*pd*line_shape;

% Check the returned spectrum against the constant-transition identity
result=test_close(result,'voitlander constant line',spec,reference,1e-10,1e-12,...
                  'subdivision of a constant transition must preserve total spherical-triangle area');
result=test_true(result,'voitlander finite real spectrum',isreal(spec)&&all(isfinite(spec(:)))&&all(spec(:)>=0),...
                 'a positive population difference and transition moment must give a finite non-negative real spectrum');

end


