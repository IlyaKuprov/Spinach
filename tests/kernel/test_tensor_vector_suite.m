% Tests tensor, vector, distribution, and relaxation utilities. Syntax:
%
%                    result=test_tensor_vector_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks Hermite splines, skew-normal density in the normal
% limit, Fokker-Planck vector reshaping helpers, correlation-function
% coefficients, tensor isotope-shift helpers, and small spin-system
% tensor extractors.
%
% ilya.kuprov@weizmann.ac.il

function result=test_tensor_vector_suite()

% State the utility target of the test
result=new_test_result('kernel/tensor_vector_suite',...
                       'Tensor, vector, and relaxation utilities',...
                       'Tensor and vector utility helpers must match exact small analytical cases.');

% Check cubic Hermite interpolation on an exactly representable parabola
spline_grid=linspace(0,1,6);
spline_obs=herm_spline(0,0,1,2,spline_grid);
spline_ref=spline_grid.^2;
result=test_close(result,'herm_spline quadratic reproduction',spline_obs,spline_ref,...
                  1e-15,1e-15,...
                  'endpoint values and derivatives for x^2 reproduce the quadratic exactly');

% Check skew normal density in the zero-skew normal-distribution limit
pdf_grid=-2:2;
pdf_obs=snormpdf(pdf_grid,1,2,0);
pdf_ref=exp(-0.5*((pdf_grid-1)/2).^2)/(2*sqrt(2*pi));
result=test_close(result,'snormpdf zero-skew normal limit',pdf_obs,pdf_ref,...
                  1e-14,1e-14,...
                  'Azzalini skew normal with alpha=0 is the ordinary normal density');

% Check phantom-to-Fokker-Planck state embedding
spin_state=[1;2];
phantom=[1 3;2 4];
fpl_obs=phan2fpl(phantom,spin_state);
fpl_ref=kron(phantom(:),spin_state);
result=test_close(result,'phan2fpl Kronecker embedding',fpl_obs,fpl_ref,...
                  1e-15,1e-15,...
                  'a spatial phantom is embedded by a Kronecker product with the spin state');

% Check painted-image extraction from a Fokker-Planck vector
paint_obs=fpl2phan(fpl_obs,[1;0],[2 2]);
result=test_close(result,'fpl2phan observable projection',paint_obs,phantom,...
                  1e-15,1e-15,...
                  'a coil selecting the first spin component recovers the original phantom');

% Check spatial averaging of a Fokker-Planck vector back to spin space
rho_obs=fpl2rho(fpl_obs,[2 2]);
rho_ref=mean(phantom(:))*spin_state;
result=test_close(result,'fpl2rho spatial average',rho_obs,rho_ref,...
                  1e-15,1e-15,...
                  'averaging over the spatial cells returns the mean phantom intensity times the spin state');

% Check isotropic rotational correlation coefficients and rates
corr_system=local_corr_system();
[weights,rates,states]=corrfun(corr_system,2,3,3,3,3);
result=test_close(result,'corrfun isotropic weight',weights{1},1/5,...
                  1e-15,1e-15,...
                  'for rank two isotropic diffusion the matching Wigner element has weight 1/(2n+1)');
result=test_close(result,'corrfun isotropic rate',rates{1},-1/corr_system.rlx.tau_c{1},...
                  1e-12,1e-12,...
                  'with rank two and second-rank tau_c, the isotropic decay rate is -1/tau_c');
result=test_true(result,'corrfun species state mask',...
                 isequal(states{1},logical([1;1;1;0])),...
                 'state membership follows non-zero basis rows in the requested chemical species');

% Check tensor isotropic-shift replacement without changing anisotropy
orig_tensor=diag([1 2 6]);
shifted=shift_iso({orig_tensor},1,5);
orig_aniso=orig_tensor-eye(3)*trace(orig_tensor)/3;
shift_aniso=shifted{1}-eye(3)*trace(shifted{1})/3;
result=test_close(result,'shift_iso isotropic part',trace(shifted{1})/3,5,...
                  1e-13,1e-13,...
                  'the requested isotropic tensor component is inserted exactly');
result=test_close(result,'shift_iso anisotropy preservation',shift_aniso,orig_aniso,...
                  1e-13,1e-13,...
                  'replacing the isotropic part preserves the rank-two anisotropy');

% Check direct coupling tensor extraction from both storage directions
spin_system=local_tensor_system();
spin_system.inter.coupling.matrix=cell(2);
spin_system.inter.coupling.matrix{1,2}=diag([1 2 3]);
spin_system.inter.coupling.matrix{2,1}=diag([4 5 6]);
result=test_close(result,'get_coupling bidirectional sum',...
                  get_coupling(spin_system,1,2),diag([5 7 9]),1e-15,1e-15,...
                  'coupling tensors stored in both spin orders are summed on extraction');

% Check g-tensor conversion from Zeeman interaction scaling
spin_system.inter.zeeman.ddscal={4*eye(3)};
spin_system.inter.gammas=2;
spin_system.tols.hbar=3;
spin_system.tols.muB=6;
result=test_close(result,'gtensorof scaling conversion',gtensorof(spin_system,1),-4*eye(3),...
                  1e-15,1e-15,...
                  'the g-tensor follows -ddscal*gamma*hbar/muB');

% Check isotropic Zeeman offset extraction in Hz
spin_system.inter.zeeman.matrix={2*pi*20*eye(3)};
spin_system.inter.basefrqs=2*pi*15;
result=test_close(result,'offsetof isotropic shift',offsetof(spin_system,1),-5,...
                  1e-14,1e-14,...
                  'the isotropic offset is the negative residual angular frequency divided by 2*pi');

end

% Minimal spin_system for corrfun
function spin_system=local_corr_system()

spin_system.bas.formalism='sphten-liouv';
spin_system.bas.basis=[1 0;0 1;1 1;0 0];
spin_system.chem.parts={1:2};
spin_system.rlx.tau_c={2e-9};

end

% Minimal spin_system for tensor extractors
function spin_system=local_tensor_system()

spin_system.comp.nspins=2;
spin_system.comp.isotopes={'1H','13C'};
spin_system.inter=struct();
spin_system.tols=struct();

end


