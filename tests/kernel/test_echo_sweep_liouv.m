% Compare echo_sweep in Hilbert and Fokker-Planck Liouville spaces.
% Syntax:
%
%                    result=test_echo_sweep_liouv()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% A zero-rank static limit tests identical physical evolution; a
% powder calculation checks that MAS changes the measured spectrum.
%
% ilya.kuprov@weizmann.ac.il

function result=test_echo_sweep_liouv()

% Announce the test target
fprintf('TESTING: Liouville echo-detected frequency sweep\n');
result=new_test_result('kernel/echo_sweep_liouv',...
                       'Liouville echo-detected frequency sweep',...
                       'The Fokker-Planck callback must preserve static Hilbert dynamics and respond to MAS.');

% P1 centre, using two representations of the same spin system
p1.orientation='111'; p1.nitrogen='14N';
[sys,inter]=diamond_p1(p1);
sys.magnet=6.9156;
bas.formalism='zeeman-hilb'; bas.approximation='none';
hilbert_system=basis(create(sys,inter),bas);
bas.formalism='zeeman-liouv';
liouville_system=basis(create(sys,inter),bas);

% Integer-step timing permits an independent static Hilbert comparison
parameters.rate=0;
parameters.axis=[1 1 1];
parameters.max_rank=0;
parameters.grid='single_crystal';
parameters.spins={'E'};
parameters.pulse_dur=40e-9;
parameters.pulse_frq=4e6;
parameters.tau=30e-9;
parameters.echo_win=40e-9;
parameters.timestep=10e-9;
parameters.sweep=2e6;
parameters.npoints=5;
parameters.verbose=0;
parameters.nphases=1;
parameters.rho0=state(hilbert_system,'Lz','E');
parameters.coil=state(hilbert_system,'L+','E');
hilbert=singlerot(hilbert_system,@echo_sweep,parameters,'esr');

% The Liouville callback has no explicit rotor-phase stack
parameters=rmfield(parameters,'nphases');
parameters.rho0=state(liouville_system,'Lz','E');
parameters.coil=state(liouville_system,'L+','E');
liouville=singlerot(liouville_system,@echo_sweep,parameters,'esr');
result=test_close(result,'static Hilbert/Liouville agreement',...
                  hilbert,liouville,1e-9*norm(hilbert),0,...
                  'echo integral and carrier offsets must agree');

% A zero interpulse delay is valid and requires no free propagation
parameters.tau=0;
parameters.nphases=1;
parameters.rho0=state(hilbert_system,'Lz','E');
parameters.coil=state(hilbert_system,'L+','E');
zero_hilbert=singlerot(hilbert_system,@echo_sweep,parameters,'esr');
parameters=rmfield(parameters,'nphases');
parameters.rho0=state(liouville_system,'Lz','E');
parameters.coil=state(liouville_system,'L+','E');
zero_liouville=singlerot(liouville_system,@echo_sweep,parameters,'esr');
result=test_close(result,'zero-delay Hilbert/Liouville agreement',...
                  zero_hilbert,zero_liouville,1e-9*norm(zero_hilbert),0,...
                  'tau=0 must skip free propagation in both formalisms');

% A small two-orientation powder exposes Fokker-Planck rotor evolution
parameters.max_rank=2;
parameters.grid='leb_1ang_rank_3';
parameters.pulse_dur=400e-9;
parameters.pulse_frq=416e3;
parameters.tau=300e-9;
parameters.echo_win=100e-9;
parameters.sweep=5e6;
static=singlerot(liouville_system,@echo_sweep,parameters,'esr');
parameters.rate=37e3;
mas=singlerot(liouville_system,@echo_sweep,parameters,'esr');
result=test_true(result,'finite static and MAS signals',...
                 all(isfinite(static(:)))&&all(isfinite(mas(:)))&&...
                 any(abs(static)>0)&&any(abs(mas)>0),...
                 'reduced P1 spectra must be nonzero and finite');
result=test_true(result,'rotor-dependent spectrum',...
                 norm(mas-static)>1e-6*norm(static),...
                 'powder MAS must change the echo spectrum');

end

