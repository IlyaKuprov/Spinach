% Tests explicit grids in both Hilbert-space rotor contexts. Syntax:
%
%                    result=test_ctx_custom_grids_hilbert()
%
% Outputs:
%
%     result - regression checks for loaded and supplied grids,
%              rotor Hamiltonian stacks, and signed weights
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_custom_grids_hilbert()

% State the context contract being tested
result=new_test_result('kernel/ctx_custom_grids_hilbert',...
                       'Custom Hilbert rotor grids',...
                       'explicit grids must preserve rotor Hamiltonian stacks and weighted signals');

% Build one anisotropic proton in Hilbert space
sys.magnet=14.1; sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]}; inter.zeeman.euler={[0.2 0.4 0.6]};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
parameters.spins={'1H'}; parameters.offset=0;
parameters.rho0=state(spin_system,'Lx','1H');
parameters.coil=parameters.rho0; parameters.decouple={};
parameters.rate=1000; parameters.axis=[1 1 1]; parameters.max_rank=1;
parameters.rate_outer=600; parameters.rate_inner=2400;
parameters.rank_outer=1; parameters.rank_inner=1;
parameters.axis_outer=[sqrt(2/3) 0 sqrt(1/3)];
parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))];
parameters.serial=true; parameters.verbose=0; parameters.sum_up=false;
parameters.grid='rep_3ang_100pts';
grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',...
                   'rep_3ang_100pts.mat'),'alphas','betas','gammas','weights');
contexts={@singlerot,@doublerot};

% Compare every rotor stack and a short density-matrix propagation
for n=1:numel(contexts)
    context=contexts{n}; parameters.grid='rep_3ang_100pts';
    [loaded,loaded_grid]=context(spin_system,@rotor_probe,parameters,'nmr');
    parameters.grid=grid;
    [supplied,supplied_grid]=context(spin_system,@rotor_probe,parameters,'nmr');
    label=func2str(context);
    result=test_true(result,[label ' returned grid'],...
                     isequaln(loaded_grid,supplied_grid),...
                     'the explicit grid is returned without modification');
    loaded_array=cell2mat(loaded); supplied_array=cell2mat(supplied);
    result=test_close(result,[label ' rotor stacks'],...
                      supplied_array(:,1:end-1),loaded_array(:,1:end-1),...
                      1e-12,1e-12,'all orientation-specific Hilbert Hamiltonians agree');
    result=test_close(result,[label ' propagated signals'],...
                      supplied_array(:,end),loaded_array(:,end),1e-12,1e-12,...
                      'short transverse density-matrix propagation agrees at every orientation');
    parameters.grid.weights=linspace(-0.25,1,numel(grid.weights));
    parameters.sum_up=true;
    weighted=context(spin_system,@rotor_probe,parameters,'nmr');
    expected=zeros(size(weighted));
    for k=1:numel(supplied)
        expected=expected+parameters.grid.weights(k)*supplied{k};
    end
    result=test_close(result,[label ' signed weights'],weighted,expected,...
                      1e-12,1e-12,'weights are used without normalisation or clipping');
    result=test_close(result,[label ' weighted signal'],weighted(end),expected(end),...
                      1e-12,1e-12,'the signed sum preserves the propagated transverse signal');
    parameters.grid.alphas=parameters.grid.alphas(:);
    parameters.grid.betas=parameters.grid.betas(:);
    parameters.grid.gammas=parameters.grid.gammas(:);
    parameters.grid.weights=parameters.grid.weights(:);
    columns=context(spin_system,@rotor_probe,parameters,'nmr');
    result=test_close(result,[label ' column vectors'],columns,weighted,...
                      1e-12,1e-12,'vector layout must not alter the quadrature');
    parameters.sum_up=false;
end

end

% Flatten the rotor Hamiltonians and propagate a transverse density matrix
function answer=rotor_probe(~,parameters,H,~,~)
stack=cell2mat(cellfun(@(x) full(x(:)),H,'UniformOutput',false));
stack=stack(:).';
rho=parameters.rho0; dt=1e-5;
for n=1:numel(H)
    P=expm(-1i*full(H{n})*dt); rho=P*rho*P';
end
answer=[stack trace(parameters.coil'*rho)];
end


