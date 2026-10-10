% Tests the double-rotor context with acquire(). Syntax:
%
%                    result=test_ctx_doublerot_acquire()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test runs a tiny anisotropic one-spin double-rotation calculation
% through doublerot() and checks the returned time-domain trace for basic
% physical and dimensional invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_doublerot_acquire()

% Announce the test target
fprintf('TESTING: Double-rotor acquire path\n');

% State the double-rotor target of the test
result=new_test_result('kernel/ctx_doublerot_acquire',...
                       'Double-rotor acquire path',...
                       'doublerot() must project states into double-rotor space and run acquire().');

% Build a one-spin anisotropic Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
bas.projections={+1};
spin_system=test_spin_system(sys,inter,bas);

% Set up a tiny double-rotation acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2000;
parameters.npoints=3;
parameters.rate_outer=800;
parameters.rate_inner=2400;
parameters.rank_outer=1;
parameters.rank_inner=1;
parameters.axis_outer=[sqrt(2/3) 0 sqrt(1/3)];
parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))];
parameters.grid='single_crystal';
parameters.serial=true;
parameters.verbose=0;

% Run the production double-rotor context
fid=doublerot(spin_system,@acquire,parameters,'nmr');

% Check the number of acquired points
result=test_close(result,'doublerot FID length',numel(fid),parameters.npoints,0,0,...
                  'acquire() should return one point for each requested time sample');

% Check that the zero-time signal survives double-rotor projection
fid_zero=parameters.coil'*parameters.rho0;
result=test_close(result,'doublerot zero-time signal',fid(1),fid_zero,1e-12,1e-12,...
                  'double-rotor projection must preserve the initial coil overlap');

% Check that the acquired trace is finite
result=test_true(result,'doublerot finite FID',all(isfinite(fid(:))),...
                 'short double-rotor propagation should not produce NaN or Inf values');

% Add finite relaxation to the representation-parity cases
inter.relaxation={'t1_t2'}; inter.equilibrium='zero'; inter.rlx_keep='secular';
inter.r1_rates={13}; inter.r2_rates={7};

% Compare full FIDs in both Liouville bases with unequal rotor ranks
for formalism={'sphten-liouv','zeeman-liouv'}
    bas.formalism=formalism{1}; bas=rmfield(bas,'projections');
    inter_form=inter;
    if strcmp(formalism{1},'zeeman-liouv')
        inter_form=rmfield(inter_form,{'r1_rates','r2_rates'});
        inter_form.relaxation={};
    end
    spin_system=test_spin_system(sys,inter_form,bas);
    parameters.rho0=state(spin_system,'L+','1H');
    parameters.coil=state(spin_system,'L+','1H');
    parameters.rank_outer=1; parameters.rank_inner=2;
    parameters.npoints=8;
    for grid={'single_crystal','rep_2ang_100pts_oct'}
        parameters.grid=grid{1};
        for rates=[800 -800 0;2400 1300 -1100]
            parameters.rate_outer=rates(1); parameters.rate_inner=rates(2);
            spin_system.sys.enable={};
            reference=doublerot(spin_system,@acquire,parameters,'nmr');
            spin_system.sys.enable={'polyadic'};
            observed=doublerot(spin_system,@acquire,parameters,'nmr');
            label=[formalism{1} '/' grid{1} '/' num2str(rates(1))];
            result=test_close(result,label,observed,reference,1e-8,1e-8,...
                              'signed rotor actions preserve the complete anisotropic FID');
        end
    end
    parameters.grid='single_crystal'; parameters.rank_outer=0;
    spin_system.sys.enable={};
    reference=doublerot(spin_system,@acquire,parameters,'nmr');
    spin_system.sys.enable={'polyadic'};
    observed=doublerot(spin_system,@acquire,parameters,'nmr');
    result=test_close(result,'zero-rank rotor',observed,reference,1e-8,1e-8,...
                      'a one-point rotor axis retains its zero derivative');
    bas.projections={+1};
end

end


