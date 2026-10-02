% Tests switch-aware periodic Fourier Laplacians. Syntax:
%
%                      result=test_fourlap_poly()
%
% Outputs:
%
%    result - Laplacian, adjoint, Nyquist, and diffusion checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_fourlap_poly()

% Describe explicit spectral Laplacian references
result=new_test_result('kernel/fourlap_poly','Implicit Fourier Laplacian',...
                       'Periodic Laplacian actions must preserve the original explicit operator.');

% Exercise unequal periods and odd/even axes in one to three dimensions
for npoints={1,8,9,[4 5],[3 4 5]}
    points=npoints{1}; periods=2+(1:numel(points));
    reference=fourlap(points,periods); nvoxels=prod(points);
    spin_system.sys.enable={};
    explicit=fourlap_poly(spin_system,points,periods);
    result=test_close(result,'disabled route',explicit,reference,0,0,...
                      'the switch-aware disabled route calls the original fourlap');
    spin_system.sys.enable={'polyadic'};
    implicit=fourlap_poly(spin_system,points,periods);
    rhs=reshape(sin(1:(2*nvoxels))+1i*cos(1:(2*nvoxels)),nvoxels,2);
    result=test_close(result,'Laplacian action',implicit*rhs,reference*rhs,...
                      1e-10,1e-11,'every axis uses its physical period and tensor ordering');
    result=test_close(result,'Laplacian adjoint',implicit'*rhs,reference'*rhs,...
                      1e-10,1e-11,'the spectral Laplacian is Hermitian');
    result=test_close(result,'constant mode',implicit*ones(nvoxels,1),zeros(nvoxels,1),...
                      1e-10,1e-11,'constants have zero periodic Laplacian');

    % Check composed dissipative propagation against an explicit exponential
    spin_system.sys.disable={}; spin_system.sys.output='hush';
    spin_system.bas.formalism='zeeman-liouv';
    spin_system.tols.prop_chop=1e-12; spin_system.tols.small_matrix=200;
    spin_system.tols.dense_matrix=0.5;
    observed=step(spin_system,1i*implicit,rhs,0.01);
    result=test_close(result,'diffusion propagation',observed,expm(0.01*full(reference))*rhs,...
                      1e-10,1e-10,'the implicit Laplacian retains diffusion decay');
end

% Preserve the even-grid second derivative of the Nyquist mode
spin_system.sys.enable={'polyadic'};
nyquist=(-1).^(0:7).';
implicit=fourlap_poly(spin_system,8,2*pi);
result=test_close(result,'Nyquist second derivative',implicit*nyquist,-16*nyquist,...
                  1e-12,1e-12,'the second derivative retains the Nyquist eigenvalue');

end


