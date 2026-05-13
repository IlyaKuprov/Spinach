% Tests compact imaging() and meshflow() context hand-off paths. Syntax:
%
%                    result=test_dynamic_fp_contexts()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that both contexts assemble finite, correctly sized
% generators and phantom-derived initial and detection states.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_fp_contexts()

% Announce the test target
fprintf('TESTING: Fokker-Planck imaging and meshflow contexts\n');

% State the context target of the test
result=new_test_result('kernel/dynamic_fp_contexts',...
                       'Fokker-Planck imaging and meshflow contexts',...
                       'imaging() and meshflow() must assemble context generators and phantom states with conserved spatial flow.');

% Exercise the Cartesian-grid imaging context
result=local_test_imaging(result);

% Exercise the unstructured-mesh flow context
result=local_test_meshflow(result);

end


function result=local_test_imaging(result)

% Build a one-spin spherical-tensor Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spn_dim=size(spin_system.bas.basis,1);

% Set a minimal one-dimensional imaging grid
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=0;
parameters.dims=0.01;
parameters.npts=10;
parameters.deriv={'period',3};
parameters.u=1e-3*ones(parameters.npts,1);
parameters.diff=1e-9;

% Supply relaxation, initial-state, and coil phantoms
parameters.rlx_ph={zeros(parameters.npts,1)};
parameters.rlx_op={sparse(spn_dim,spn_dim)};
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H')};

% Run the production context into a probe pulse sequence
answer=imaging(spin_system,@local_context_probe,parameters);
problem_dim=parameters.npts*spn_dim;

% Check dimensions of all assembled generators
result=test_close(result,'imaging spin dimension',answer.spn_dim,spn_dim,0,0,...
                  'the context must pass the basis dimension to the pulse sequence');
result=test_close(result,'imaging problem dimension',answer.problem_dim,problem_dim,0,0,...
                  'the Cartesian spatial dimension times spin dimension gives the Fokker-Planck size');
result=local_test_square(result,'imaging H size',answer.H,problem_dim);
result=local_test_square(result,'imaging R size',answer.R,problem_dim);
result=local_test_square(result,'imaging K size',answer.K,problem_dim);
result=local_test_square(result,'imaging Fx size',answer.F,problem_dim);
result=local_test_square(result,'imaging Gx size',answer.G{1},problem_dim);
result=test_true(result,'imaging transverse gradients',isempty(answer.G{2})&&isempty(answer.G{3}),...
                 'one-dimensional imaging grids must leave transverse gradient operators empty');

% Check phantom state construction and flow conservation
result=test_close(result,'imaging rho0 length',numel(answer.rho0),problem_dim,0,0,...
                  'one spin state must be placed in every imaging voxel');
result=test_close(result,'imaging coil length',numel(answer.coil),problem_dim,0,0,...
                  'one detection state must be placed in every imaging voxel');
result=test_close(result,'imaging flow conservation',sum(full(answer.F),1),zeros(1,problem_dim),1e-10,1e-10,...
                  'periodic one-dimensional flow and diffusion must conserve total spatial mass');

end


function result=local_test_meshflow(result)

% Build a one-spin spherical-tensor Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spn_dim=size(spin_system.bas.basis,1);

% Attach a two-cell Voronoi mesh with one shared boundary
spin_system.mesh=local_two_cell_mesh();
spc_dim=spin_system.mesh.vor.ncells;
problem_dim=spc_dim*spn_dim;

% Supply operator phantoms and state phantoms
parameters.H_ph={ones(spc_dim,1)};
parameters.H_op={sparse(spn_dim,spn_dim)};
parameters.R_ph={zeros(spc_dim,1)};
parameters.R_op={sparse(spn_dim,spn_dim)};
parameters.K_ph={zeros(spc_dim,1)};
parameters.K_op={speye(spn_dim)};
parameters.rho0_ph={ones(spc_dim,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(spc_dim,1)};
parameters.coil_st={state(spin_system,'Lz','1H')};
parameters.diff=1e-8;

% Run the production meshflow context into a probe pulse sequence
answer=meshflow(spin_system,@local_context_probe,parameters);

% Check dimensions of assembled unstructured-mesh generators
result=test_close(result,'meshflow spin dimension',answer.spn_dim,spn_dim,0,0,...
                  'meshflow() must pass the basis dimension to the pulse sequence');
result=test_close(result,'meshflow problem dimension',answer.problem_dim,problem_dim,0,0,...
                  'the mesh cell count times spin dimension gives the Fokker-Planck size');
result=local_test_square(result,'meshflow H size',answer.H,problem_dim);
result=local_test_square(result,'meshflow R size',answer.R,problem_dim);
result=local_test_square(result,'meshflow K size',answer.K,problem_dim);
result=local_test_square(result,'meshflow F size',answer.F,problem_dim);
result=test_close(result,'meshflow rho0 length',numel(answer.rho0),problem_dim,0,0,...
                  'one spin state must be placed in every active Voronoi cell');
result=test_close(result,'meshflow coil length',numel(answer.coil),problem_dim,0,0,...
                  'one detection state must be placed in every active Voronoi cell');
result=test_close(result,'meshflow flow conservation',sum(full(answer.F),1),zeros(1,problem_dim),1e-15,1e-15,...
                  'finite-volume mesh diffusion with closed boundaries must conserve total mass');

end


function answer=local_context_probe(~,parameters,H,R,K,G,F)

% Return context products for regression checks
answer.spc_dim=parameters.spc_dim;
answer.spn_dim=parameters.spn_dim;
answer.problem_dim=parameters.spc_dim*parameters.spn_dim;
answer.rho0=parameters.rho0;
answer.coil=parameters.coil;
answer.H=H;
answer.R=R;
answer.K=K;
answer.G=G;
answer.F=F;

end


function result=local_test_square(result,label,A,matrix_dim)

% Check that a context object is a square matrix of the expected size
result=test_close(result,label,size(A),[matrix_dim matrix_dim],0,0,...
                  'context generators must act on the full Fokker-Planck state vector');
result=test_true(result,[label ' finite'],all(isfinite(full(A(:)))),...
                 'assembled context generators must not contain NaN or Inf values');

end


function mesh=local_two_cell_mesh()

% Active mesh vertices are the parents of the two finite-volume cells
mesh.idx.active=[1 2];
mesh.idx.triangles=[1 2 3];

% Voronoi cells share the boundary between vertices one and two
mesh.vor.ncells=2;
mesh.vor.weights=[1;1];
mesh.vor.cells={[1 2 3],[2 1 4]};
mesh.vor.vertices=[0 0;0 1;-1 0;1 0];

% Cell-centre coordinates and zero velocity field
mesh.x=[0;1];
mesh.y=[0;0];
mesh.u=[0;0];
mesh.v=[0;0];

end


