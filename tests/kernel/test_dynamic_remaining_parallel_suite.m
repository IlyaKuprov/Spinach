% Tests remaining parallel, stochastic, and diagnostic utilities. Syntax:
%
%                    result=test_dynamic_remaining_parallel_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks distributed-array reconstruction, zero stochastic
% Redfield integration, and Fokker-Planck overwinding diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_remaining_parallel_suite()

% State the utility target of the test
result=new_test_result('kernel/dynamic_remaining_parallel_suite',...
                       'Remaining parallel and stochastic utilities',...
                       'Parallel and stochastic helper utilities must preserve compact reference cases.');

% Keep compact parallel smoke paths to one local worker on this host
local_ensure_pool();

% Check dimension-specific distributed array construction when the toolbox is present
if (exist('distributed','class')==8)&&(exist('codistributor1d','class')==8)
    dense_array=reshape(1:12,[3 4]);
    distributed_array=distrib_dim(dense_array,2);
    result=test_close(result,'distrib_dim reconstruction',gather(distributed_array),dense_array,...
                      1e-14,1e-14,...
                      'gathering a codistributed array must reproduce the original dense array');
else
    result=test_true(result,'distrib_dim toolbox unavailable',true,...
                     'distributed-array coverage is blocked because the Parallel Computing Toolbox is unavailable');
end

% Check numerical Redfield integration gives zero relaxation for zero stochastic Hamiltonians
spin_system=local_liouvillian_system(1);
H0=sparse(1,1,1e-3,1,1);
H1=repmat({sparse(1,1)},2001,1);
R=ngce(spin_system,H0,H1,1,10,0);
result=test_close(result,'ngce zero stochastic Hamiltonian',R,sparse(1,1),1e-14,1e-14,...
                  'a zero stochastic Hamiltonian trajectory must integrate to a zero relaxation superoperator');

% Check overwinding diagnostics complete and draw a spectrum for a safe one-dimensional grid
figures_before=numel(findall(0,'Type','figure'));
rho=ones(10,1);
overwound(rho,[10 1 1],1);
figures_after=numel(findall(0,'Type','figure'));
close all;
result=test_true(result,'overwound one-dimensional diagnostic',figures_after>=figures_before+1,...
                 'a resolved one-dimensional Fokker-Planck state should produce one diagnostic spectrum');

end


function spin_system=local_liouvillian_system(dim)

% Create a quiet spherical-tensor Liouville descriptor
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.bas.formalism='sphten-liouv';
spin_system.bas.basis=zeros(dim,1);
spin_system.tols.liouv_zero=1e-12;
spin_system.tols.prop_chop=1e-12;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.small_matrix=10;

end


function local_ensure_pool()

% Start a one-worker process pool for compact parallel utilities
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',1);
end

end

