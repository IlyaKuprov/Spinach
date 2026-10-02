% Tests implicit Fourier steady-state DNP scans. Syntax:
%
%                         result=test_dnp_fft()
%
% Outputs:
%
%    result - steady-state observables and explicit-route checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_dnp_fft()

% State the independent explicit steady-state reference
result=new_test_result('kernel/dnp_fft','FFT DNP steady states',...
                       'Implicit GMRES scans must reproduce explicit backslash observables.');

% Build an electron-proton pair with a finite unthermalised decay generator
sys.magnet=0.0001; sys.isotopes={'E','1H'};
inter.zeeman.scalar={2.0023,0};
inter.coupling.scalar=cell(2); inter.coupling.scalar{1,2}=1e4;
bas.approximation='none';

% Exercise both Liouville bases, signed offsets, phase counts, and coils
for formalism={'sphten-liouv','zeeman-liouv'}
    bas.formalism=formalism{1};
    spin_system=test_spin_system(sys,inter,bas);
    spin_system=assume(spin_system,'labframe');
    H=hamiltonian(spin_system); R=-1000*speye(size(H)); K=kinetics(spin_system);
    parameters.mw_pwr=2*pi*2e4; parameters.mw_frq=2*pi*[-1e4 0 1e4];
    parameters.g_ref=2.0023;
    parameters.rho0=state(spin_system,'Lz','E');
    parameters.coil=[state(spin_system,'Lz','E') state(spin_system,'Lz','1H')];
    parameters.mw_oper=operator(spin_system,'Lx','E');
    for nphases=[1 4 5 8]
        parameters.nphases=nphases;
        spin_system.sys.enable={}; parameters.method='fp-backs';
        reference=dnp_freq_scan(spin_system,parameters,H,R,K);
        parameters.method='fp-gmres';
        spin_system.sys.enable={'polyadic'};
        observed=dnp_freq_scan(spin_system,parameters,H,R,K);
        label=[formalism{1} '/' num2str(nphases)];
        result=test_close(result,['implicit ' label],observed,reference,1e-8,1e-8,...
                          'complete multi-coil frequency scans match direct solves');
        parameters.method='fp-backs';
        unchanged=dnp_freq_scan(spin_system,parameters,H,R,K);
        result=test_close(result,['backslash ' label],unchanged,reference,0,0,...
                          'direct FP solves remain explicit even when polyadics are enabled');
    end

    % Check the original matrix GMRES against an exact zero-drive limit
    zero_params=parameters; zero_params.nphases=1; zero_params.mw_pwr=0;
    zero_params.method='fp-gmres'; spin_system.sys.enable={};
    observed=dnp_freq_scan(spin_system,zero_params,0*H,R,0*K);
    reference=repmat((parameters.coil'*parameters.rho0).',numel(parameters.mw_frq),1);
    result=test_close(result,'explicit GMRES retained',observed,reference,1e-10,1e-10,...
                      'zero microwave drive retains the supplied equilibrium polarisation');

    % Compare the new iterative route at a larger physical carrier frequency
    sys_high=sys; sys_high.magnet=0.005;
    high_system=test_spin_system(sys_high,inter,bas);
    high_system=assume(high_system,'labframe'); H_high=hamiltonian(high_system);
    R_high=-1000*speye(size(H_high)); K_high=kinetics(high_system);
    high_params=parameters; high_params.mw_frq=2*pi*[-1e6 0 1e6];
    high_params.method='fp-backs'; high_system.sys.enable={};
    reference=dnp_freq_scan(high_system,high_params,H_high,R_high,K_high);
    high_params.method='fp-gmres'; high_system.sys.enable={'polyadic'};
    observed=dnp_freq_scan(high_system,high_params,H_high,R_high,K_high);
    result=test_close(result,'larger carrier frequency',observed,reference,1e-8,1e-8,...
                      'Fourier-block preconditioning retains physical high-frequency steady states');

    % Check that the LvN direct solver is not changed by the enable switch
    spin_system=assume(spin_system,'esr'); H=hamiltonian(spin_system);
    parameters.ez_oper=operator(spin_system,'Lz','E'); parameters.method='lvn-backs';
    spin_system.sys.enable={}; reference=dnp_freq_scan(spin_system,parameters,H,R,K);
    spin_system.sys.enable={'polyadic'}; observed=dnp_freq_scan(spin_system,parameters,H,R,K);
    result=test_close(result,'LvN unchanged',observed,reference,0,0,...
                      'the rotating-frame direct solver is unaffected');
end

end


