% Compare small MAS echoes with explicit and FFT rotor derivatives.
% Syntax:
%
%                         result=test_fft_echo()
%
% Outputs:
%
%    result - complex echo parity checks in both Liouville formalisms
%
% ilya.kuprov@weizmann.ac.il

function result=test_fft_echo()

% State the comparison without changing the rotor rank
result=new_test_result('kernel/fft_echo','FFT rotor echoes',...
                       'Implicit and explicit same-rank rotor generators must give identical complex echoes.');

% A small anisotropic two-spin system exercises noncommuting MAS dynamics
sys.isotopes={'1H','13C'}; sys.magnet=1;
sys.parallel={'processes',2};
inter.zeeman.scalar={1 2}; inter.zeeman.eigs={[0 10 20],[0 0 0]};
inter.zeeman.euler={[0 0 0],[0 0 0]};
inter.coupling.scalar={0 100;100 0};
bas.approximation={'none'};
parameters.axis=[1 1 1]; parameters.grid='leb_2ang_rank_5';
parameters.spins={'1H'}; parameters.offset=0;
parameters.pulse_dur=1e-4; parameters.pulse_frq=2500;
parameters.tau=2e-4; parameters.echo_win=1e-4;
parameters.timestep=2e-5; parameters.sweep=1000; parameters.npoints=3;
parameters.verbose=0; parameters.serial=true;

% Both spin bases use the same static and signed spinning specifications
for formalism={'zeeman-liouv','sphten-liouv'}
    bas.formalism=formalism{1};
    spin_system=basis(create(sys,inter),bas);
    parameters.rho0=state(spin_system,'Lz','1H');
    parameters.coil=state(spin_system,'L+','1H');
    for rotor_case=[0 1 2]
        parameters.max_rank=rotor_case;
        parameters.serial=(rotor_case~=2);
        parameters.rate=(-1)^rotor_case*500;
        spin_system.sys.enable={};
        reference=singlerot(spin_system,@echo_sweep,parameters,'nmr');
        spin_system.sys.enable={'polyadic'};
        implicit=singlerot(spin_system,@echo_sweep,parameters,'nmr');
        result=test_close(result,[formalism{1} ' rank ' num2str(rotor_case)],...
                          implicit,reference,1e-10*norm(reference),0,...
                          'complete complex echoes agree at identical rank and carrier offsets');
    end
end


% Analytical decoupling remains available for implicit acquisition
acq_par=parameters; acq_par.grid='single_crystal';
acq_par.max_rank=1; acq_par.rate=500; acq_par.serial=true;
acq_par.npoints=4; acq_par.decouple={'13C'};
acq_par.rho0=state(spin_system,'L+','1H'); acq_par.coil=acq_par.rho0;
spin_system.sys.enable={};
reference=singlerot(spin_system,@acquire,acq_par,'nmr');
spin_system.sys.enable={'polyadic'};
implicit=singlerot(spin_system,@acquire,acq_par,'nmr');
result=test_close(result,'decoupled acquisition',implicit,reference,...
                  1e-10*norm(reference),0,...
                  'projected implicit propagation matches analytical decoupling');

% GPU acquisition exercises FFT actions and analytical decoupling
if gpuDeviceCount('available')>0
    spin_system.sys.enable={'gpu'};
    gpu_exp=singlerot(spin_system,@acquire,acq_par,'nmr');
    spin_system.sys.enable={'gpu','polyadic'};
    gpu_fft=singlerot(spin_system,@acquire,acq_par,'nmr');
    result=test_close(result,'GPU explicit acquisition',gather(gpu_exp),reference,...
                      1e-10*norm(reference),0,'GPU propagation agrees with the CPU reference');
    result=test_close(result,'GPU FFT acquisition',gather(gpu_fft),reference,...
                      1e-10*norm(reference),0,'FFT propagation preserves decoupling on the GPU');
else
    result.messages{end+1}='SKIP: GPU acquisitions require a usable GPU.';
end

% The motivating P1 ESR model uses the same small one-orientation grid
p1.orientation='111'; p1.nitrogen='14N';
[sys,inter]=diamond_p1(p1);
sys.magnet=6.9156; sys.parallel={'processes',2};
bas.formalism='zeeman-liouv';
spin_system=basis(create(sys,inter),bas);
parameters.grid='single_crystal'; parameters.max_rank=2; parameters.serial=true;
parameters.spins={'E'}; parameters.rate=37e3;
parameters.pulse_dur=400e-9; parameters.pulse_frq=416e3;
parameters.tau=300e-9; parameters.echo_win=100e-9;
parameters.timestep=10e-9; parameters.sweep=5e6;
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
reference=singlerot(spin_system,@echo_sweep,parameters,'esr');
spin_system.sys.enable={'polyadic'};
implicit=singlerot(spin_system,@echo_sweep,parameters,'esr');
result=test_close(result,'P1 ESR echo',implicit,reference,...
                  1e-10*norm(reference),0,...
                  'finite-pulse P1 complex echoes agree at identical rotor rank');

end


