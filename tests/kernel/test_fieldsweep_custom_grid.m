% Tests named and explicit initial grids through adaptive fieldsweep. Syntax:
%
%                    result=test_fieldsweep_custom_grid()
%
% Outputs:
%
%     result - regression checks for anisotropic spectra, field axes,
%              and adaptive spherical-triangle subdivision
%
% ilya.kuprov@weizmann.ac.il

function result=test_fieldsweep_custom_grid()

% State the shared adaptive integration contract
result=new_test_result('kernel/fieldsweep_custom_grid',...
                       'Custom-grid adaptive field spectra',...
                       'An explicit initial grid must use the named-grid adaptive algorithm.');

% Build a rhombic electron whose transition fields vary across the sphere
sys.magnet=1; sys.isotopes={'E'};
inter.zeeman.eigs={[2 2.1 2.2]}; inter.zeeman.euler={[0 0 0]};
bas.formalism='zeeman-hilb'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
parameters.spins={'E'}; parameters.mw_freq=9.5e9;
parameters.window=[0.28 0.38]; parameters.fwhm=0.005;
parameters.npoints=31; parameters.tm_tol=1e-10;
parameters.rspt_order=Inf; parameters.int_tol=Inf;
grid_name='icos_2ang_12pts';
sph_grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',grid_name));

% Compare identical initial grids without recursive refinement
parameters.grid=grid_name;
[coarse_named,namedpar]=fieldsweep(spin_system,parameters);
parameters.grid=sph_grid;
[coarse_struct,structpar]=fieldsweep(spin_system,parameters);
result=test_close(result,'initial-grid spectrum',coarse_struct,coarse_named,...
                  1e-12,1e-12,'the same vertices enter the same triangular integration');
result=test_close(result,'initial-grid field axis',structpar.b_axis,namedpar.b_axis,...
                  0,0,'the grid representation does not change the field axis');

% Require recursive refinement and compare both grid representations again
parameters.int_tol=1e-5; parameters.grid=grid_name;
[refined_named,namedpar]=fieldsweep(spin_system,parameters);
parameters.grid=sph_grid;
[refined_struct,structpar]=fieldsweep(spin_system,parameters);
result=test_close(result,'adaptive spectrum',refined_struct,refined_named,...
                  1e-12,1e-12,'explicit grids retain adaptive spherical-triangle subdivision');
result=test_close(result,'adaptive field axis',structpar.b_axis,namedpar.b_axis,...
                  0,0,'adaptive subdivision preserves the common field axis');
result=test_true(result,'nonzero anisotropic spectrum',norm(refined_named)>0,...
                 'the rhombic electron produces a nonzero field-swept signal');
result=test_true(result,'adaptive subdivision used',...
                 norm(refined_named-coarse_named)>sqrt(eps)*norm(refined_named),...
                 'tightening the recursion criterion changes the integrated spectrum beyond roundoff');
fprintf('FIELDSWEEP_ADAPTIVE_DELTA %.17g NAMED_STRUCT_ERROR %.17g\n',...
        norm(refined_named-coarse_named),norm(refined_struct-refined_named));

end


