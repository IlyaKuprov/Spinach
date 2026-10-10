% Tests implicit Fourier flow and diffusion operators. Syntax:
%
%                       result=test_hydro_fft()
%
% Outputs:
%
%    result - spatial actions, adjoints, and diffusion checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_hydro_fft()

% Define independently materialised reference actions
result=new_test_result('kernel/hydro_fft','FFT hydrodynamics',...
                       'Implicit spatial derivatives must preserve flow and diffusion actions.');

% Set the spin-space lift and numerical cleanup tolerance
spin_system.bas.basis={speye(2)}; spin_system.bas.offsets=[0;2]; spin_system.tols.liouv_zero=1e-14;
spin_system.tols.dense_matrix=0.5; spin_system.tols.small_matrix=200; spin_system.sys.disable={};

% Exercise tensor ordering and odd/even grids in one to three dimensions
for npts={11,[10 11],[10 11 12]}
    parameters=struct(); parameters.npts=npts{1}; parameters.dims=1:numel(parameters.npts);
    parameters.deriv={'fourier'};
    nvoxels=prod(parameters.npts);
    rhs=reshape(sin(1:(2*nvoxels))+1i*cos(1:(2*nvoxels)),nvoxels,2);
    spin_system.sys.enable={};
    [ref_x,ref_y,ref_z]=hydrodynamics(spin_system,parameters);
    spin_system.sys.enable={'polyadic'};
    [obs_x,obs_y,obs_z]=hydrodynamics(spin_system,parameters);
    refs={ref_x,ref_y,ref_z}; observed={obs_x,obs_y,obs_z};
    for axis=1:numel(parameters.npts)
        label=[num2str(numel(parameters.npts)) '/' num2str(axis)];
        result=test_close(result,['action ' label],observed{axis}*rhs,refs{axis}*rhs,...
                          1e-10,1e-11,'the implicit axis action preserves tensor ordering');
        result=test_close(result,['adjoint ' label],observed{axis}'*rhs,refs{axis}'*rhs,...
                          1e-10,1e-11,'the Hermitian momentum generator retains its adjoint');
        result=test_close(result,['constant ' label],observed{axis}*ones(nvoxels,1),...
                          zeros(nvoxels,1),1e-10,1e-11,'a periodic constant has zero derivative');
    end
    result=test_true(result,'absent axes',all(cellfun(@isempty,observed(numel(parameters.npts)+1:end))),...
                     'absent spatial dimensions keep empty outputs');

    % Compare composed variable-coefficient flow and anisotropic diffusion
    grid_shape=parameters.npts;
    if isscalar(grid_shape), grid_shape=[grid_shape 1]; end %#ok<AGROW>
    parameters.u=reshape(0.02+0.01*sin(1:nvoxels),grid_shape);
    parameters.diff=0.1*eye(numel(parameters.npts))+0.01*ones(numel(parameters.npts),numel(parameters.npts));
    spin_system.sys.enable={}; reference=v2fplanck(spin_system,parameters);
    spin_system.sys.enable={'polyadic'}; implicit=v2fplanck(spin_system,parameters);
    states=[rhs;rhs];
    result=test_close(result,'flow/diffusion action',implicit*states,reference*states,...
                      1e-9,1e-10,'variable coefficients and first-derivative products are preserved');
    result=test_close(result,'flow/diffusion adjoint',implicit'*states,reference'*states,...
                      1e-9,1e-10,'composed spatial dynamics retain their adjoint');

    % Exercise nonuniform diffusion independently of the uniform tensor
    parameters=rmfield(parameters,'diff');
    parameters.dxx=reshape(0.01+0.001*cos(1:nvoxels),grid_shape);
    if numel(parameters.npts)>1
        parameters.dxy=zeros(grid_shape); parameters.dyx=parameters.dxy;
        parameters.dyy=parameters.dxx;
    end
    if numel(parameters.npts)>2
        parameters.dxz=zeros(grid_shape); parameters.dyz=parameters.dxz;
        parameters.dzx=parameters.dxz; parameters.dzy=parameters.dxz;
        parameters.dzz=parameters.dxx;
    end
    spin_system.sys.enable={}; reference=v2fplanck(spin_system,parameters);
    spin_system.sys.enable={'polyadic'}; implicit=v2fplanck(spin_system,parameters);
    result=test_close(result,'nonuniform diffusion',implicit*states,reference*states,...
                      1e-9,1e-10,'derivative-coefficient-derivative products retain their action');
end

% Confirm that finite-difference representations are unchanged
parameters=struct('npts',10,'dims',1,'deriv',{{'period',5}});
spin_system.sys.enable={}; [reference,~,~]=hydrodynamics(spin_system,parameters);
spin_system.sys.enable={'polyadic'}; [implicit,~,~]=hydrodynamics(spin_system,parameters);
result=test_close(result,'finite difference unchanged',inflate(implicit),reference,...
                  0,0,'periodic finite differences still have numeric polyadic cores');

end


