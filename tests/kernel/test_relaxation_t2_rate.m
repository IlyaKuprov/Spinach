% Tests phenomenological T2 relaxation rate. Syntax:
%
%                    result=test_relaxation_t2_rate()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that the t1_t2 relaxation model gives transverse L+ order
% the negative generator eigenvalue corresponding to the specified R2 rate.
%
% ilya.kuprov@weizmann.ac.il

function result=test_relaxation_t2_rate()

% State the relaxation target of the test
result=new_test_result('kernel/relaxation_t2_rate',...
                       'Phenomenological T2 decay rate',...
                       'the t1_t2 model must assign transverse magnetisation the generator eigenvalue -R2.');

% Build a one-spin relaxation system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
inter.relaxation={'t1_t2'};
inter.r1_rates={2.0};
inter.r2_rates={7.0};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.temperature=298;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Build relaxation superoperator and transverse state
R=relaxation(spin_system);
rho=state(spin_system,'L+','1H');
Rrho=R*rho;

% Check that L+ is an eigenstate with the negative R2 generator eigenvalue
result=test_close(result,'R2 eigenvalue on L+',Rrho,-7.0*rho,1e-12,1e-12,...
                  'decay generators have negative eigenvalues, with magnitude R2');

end

