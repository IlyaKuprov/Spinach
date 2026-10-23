% Tests direct-grid field spectra against analytic Lorentzian lines. Syntax:
%
%                    result=test_fieldsweep_custom_grid()
%
% Outputs:
%
%     result - regression checks for line shape, normalisation,
%              output dimensions, and signed orientation weights
%
% ilya.kuprov@weizmann.ac.il

function result=test_fieldsweep_custom_grid()

% State the direct quadrature contract being tested
result=new_test_result('kernel/fieldsweep_custom_grid',...
                       'Custom-grid field spectra',...
                       'direct quadrature must preserve analytic single-spin lines and signed sums');

% Build an electron with a single analytically known transition
sys.magnet=1; sys.isotopes={'E'};
inter.zeeman.eigs={[2 2.1 2.2]}; inter.zeeman.euler={[0 0 0]};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
parameters.spins={'E'}; parameters.mw_freq=9.5e9;
parameters.window=[0.28 0.38]; parameters.fwhm=0.002;
parameters.npoints=101; parameters.tm_tol=1e-10;
parameters.rspt_order=Inf; parameters.int_tol=1e-3;
parameters.grid.alphas=0; parameters.grid.betas=0;
parameters.grid.gammas=0; parameters.grid.weights=1.75;
parameters.rho0=[0 0;0 1];

% Evaluate the orientation along the principal z axis
[observed,returned]=fieldsweep(spin_system,parameters);
[I,Q]=hamiltonian(assume(spin_system,'labframe','zeeman'));
Hz=full(I+orientation(Q,[0 0 0]));
slope=abs(Hz(1,1)-Hz(2,2));
centre=2*pi*parameters.mw_freq/slope;
Hmw=full(operator(spin_system,'Lx','E'));
intensity=parameters.grid.weights*abs(Hmw(1,2))^2*abs(spin('E'))/slope;
gamma=parameters.fwhm/2;
expected=(intensity/(pi*gamma))./...
         (1+((returned.b_axis-centre)/gamma).^2);
result=test_close(result,'analytic single orientation',observed,expected,...
                  1e-9,1e-9,'the line has the analytic field centre and unit-area Lorentzian amplitude');
result=test_true(result,'row spectrum dimensions',...
                 isequal(size(observed),[1 parameters.npoints]),...
                 'orientation accumulation returns one spectrum row, not an outer product');
result=test_true(result,'field axis dimensions',...
                 isequal(size(returned.b_axis),size(observed)),...
                 'the returned axis and spectrum have identical dimensions');

% Compare an anisotropic signed sum with independent one-point calls
parameters.grid.alphas=[0 0 0];
parameters.grid.betas=[0 pi/2 pi/2];
parameters.grid.gammas=[0 0 pi/2];
parameters.grid.weights=[1.75 -0.6 0.2];
combined=fieldsweep(spin_system,parameters);
expected=zeros(size(combined));
for n=1:numel(parameters.grid.weights)
    local=parameters;
    local.grid.alphas=parameters.grid.alphas(n);
    local.grid.betas=parameters.grid.betas(n);
    local.grid.gammas=parameters.grid.gammas(n);
    local.grid.weights=parameters.grid.weights(n);
    expected=expected+fieldsweep(spin_system,local);
end
result=test_close(result,'signed orientation sum',combined,expected,...
                  1e-10,1e-10,'negative quadrature weights are retained in a linear sum');

% Check column layouts and an overall negative scale factor
parameters.grid.alphas=parameters.grid.alphas(:);
parameters.grid.betas=parameters.grid.betas(:);
parameters.grid.gammas=parameters.grid.gammas(:);
parameters.grid.weights=-2*parameters.grid.weights(:);
observed=fieldsweep(spin_system,parameters);
result=test_close(result,'negative scale and columns',observed,-2*combined,...
                  1e-10,1e-10,'weights are neither normalised nor clipped');

end


