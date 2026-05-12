% Tests powder averaging against explicit weighted summation. Syntax:
%
%                    result=test_ctx_powder_average()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test asks powder() for individual orientation traces and checks that
% the default powder average is the same weighted sum.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_powder_average()

% State the powder-averaging target of the test
result=new_test_result('kernel/ctx_powder_average',...
                       'Powder weighted sum path',...
                       'powder() must sum orientation outputs with grid weights.');

% Build a one-spin anisotropic Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
spin_system=test_spin_system(sys,inter,bas);

% Set up a tiny powder acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2000;
parameters.npoints=3;
parameters.grid='leb_2ang_rank_5';
parameters.serial=true;
parameters.verbose=0;

% Run the averaged powder calculation
[fid_avg,sph_grid]=powder(spin_system,@acquire,parameters,'nmr');

% Run the per-orientation powder calculation
parameters.sum_up=false;
fid_cells=powder(spin_system,@acquire,parameters,'nmr');

% Assemble the independent weighted sum
fid_ref=sph_grid.weights(1)*fid_cells{1};
for n=2:numel(fid_cells)
    fid_ref=fid_ref+sph_grid.weights(n)*fid_cells{n};
end

% Check the powder average accumulation path
result=test_close(result,'weighted powder sum',fid_avg,fid_ref,1e-12,1e-12,...
                  'the powder average must equal the explicit grid-weighted sum of all orientations');

end


