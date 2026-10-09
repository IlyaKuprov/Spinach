% Slow passage spectrum with several substances: every substance has
% its own stationary unit state, and the frequency grid contains the
% exact zero frequency where the unit states used to make the linear
% solve singular. The slowpass spectrum must be finite on the whole
% grid and must agree with the FFT of the acquired FID at the peaks.
% Also checks spatial embedding, selective reaction coupling, and the
% unchanged wavefunction resolvent. Syntax:
%
%                    result=test_slowpass_unit_states()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=test_slowpass_unit_states.m>

function result=test_slowpass_unit_states()

% Test record
result=new_test_result('kernel/slowpass_unit_states',...
                       'Slowpass unit states per substance',...
                       'slowpass() must stay finite and FFT-consistent when several stationary unit states are present.');

% Two one-spin substances, no reactions, zero equilibrium
sys.magnet=14.1; sys.isotopes={'1H','1H'};
inter.chem.parts={1,2}; inter.chem.concs=[1 1];
inter.zeeman.scalar={0,0.05}; inter.coupling.scalar=cell(2);
inter.relaxation={'damp'}; inter.damp_rate=8.0;
inter.equilibrium='zero'; inter.rlx_keep='labframe';
inter.temperature=298;
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system); R=relaxation(spin_system);
K=kinetics(spin_system);

% Acquisition with a grid that contains zero frequency
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={}; parameters.sweep=4096;
parameters.npoints=4096;
fid=acquire(spin_system,parameters,H,R,K);
spectrum_fft=fftshift(fft(fid));
frq_axis=ft_axis(0,parameters.sweep,parameters.npoints);
parameters.sweep=[frq_axis(1) frq_axis(end)];
spectrum_slow=slowpass(spin_system,parameters,H,R,K);

% Every grid point must be finite
result=test_true(result,'slowpass finite on a grid containing zero frequency',...
                 all(isfinite(spectrum_slow)),...
                 'stationary unit states must not make the frequency domain solve singular');

% Compare the peaks allowing for the discrete-transform baseline
zero_idx=parameters.npoints/2+1;
second_frq=1e-6*inter.zeeman.scalar{2}*spin_system.inter.basefrqs(2)/(2*pi);
[~,second_idx]=min(abs(frq_axis-second_frq)); baseline=abs(fid(1))/2;
result=test_close(result,'slowpass FFT amplitude at zero frequency',...
                  spectrum_slow(zero_idx),spectrum_fft(zero_idx),2*baseline,2e-3,...
                  'the unit state projection must leave the traceless spectrum unchanged');
result=test_close(result,'slowpass FFT amplitude at the second resonance',...
                  spectrum_slow(second_idx),spectrum_fft(second_idx),2*baseline,2e-3,...
                  'the second substance resonance must match the FFT of the same damped FID');

% Check a damped anisotropic spin in the gridfree spatial basis
clear sys inter bas parameters;
sys.magnet=14.1; sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]}; inter.zeeman.euler={[0 0 0]};
inter.relaxation={'damp'}; inter.damp_rate=8;
inter.equilibrium='zero'; inter.rlx_keep='labframe';
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
parameters.spins={'1H'}; parameters.offset=0;
parameters.rho0=state(spin_system,'L+','1H'); parameters.coil=parameters.rho0;
parameters.decouple={}; parameters.sweep=[-100 100]; parameters.npoints=3;
parameters.max_rank=2; parameters.tau_c=1e-3; parameters.verbose=0;
spectrum_sle=gridfree(spin_system,@slowpass,parameters,'nmr');
result=test_true(result,'gridfree slowpass finite',all(isfinite(spectrum_sle)),...
                 'spatially expanded unit directions must match the gridfree Liouvillian');

% Compare selective singlet and triplet loss against the unchanged resolvent
clear sys inter bas parameters;
sys.magnet=1; sys.isotopes={'E','E'};
inter.chem.parts={1:2}; inter.chem.concs=1;
inter.chem.reactions={struct('reactants',1,'products',[],'matching',zeros(0,2),...
    'rate',2,'selector',{{'singlet',[1 2]}}),...
    struct('reactants',1,'products',[],'matching',zeros(0,2),...
    'rate',3,'selector',{{'triplet',[1 2]}})};
inter.relaxation={'damp'}; inter.damp_rate=8;
inter.equilibrium='zero'; inter.rlx_keep='labframe';
bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
K=kinetics(spin_system); R=relaxation(spin_system); H=sparse(size(K,1),size(K,2));
parameters.rho0=state(spin_system,{'Lz','Lz'},{1,2});
parameters.coil=parameters.rho0; parameters.sweep=[-1 1]; parameters.npoints=3;
spectrum_rx=slowpass(spin_system,parameters,H,R,K);
reference=zeros(3,1); freq_grid=2*pi*linspace(-1,1,3);
for n=1:3
    reference(n)=3*parameters.coil'*((-R-K+1i*freq_grid(n)*speye(size(K)))\parameters.rho0);
end
result=test_close(result,'selective reaction resolvent',spectrum_rx,reference,1e-12,1e-12,...
                  'identity and spin order must retain their physical reaction coupling');

% Preserve the direct wavefunction resolvent without calling unit_state
clear sys bas parameters;
sys.magnet=1; sys.isotopes={'1H'};
bas.formalism='zeeman-wavef'; bas.approximation={'none'};
spin_system=test_spin_system(sys,struct(),bas);
parameters.rho0=[1;0]; parameters.coil=parameters.rho0;
parameters.sweep=[-1 1]; parameters.npoints=3;
spectrum_wf=slowpass(spin_system,parameters,sparse(diag([1 -1])),-speye(2),sparse(2,2));
reference=3./(1+1i*(1+2*pi*linspace(-1,1,3)'));
result=test_close(result,'wavefunction resolvent',spectrum_wf,reference,1e-12,1e-12,...
                  'wavefunction inputs have no Liouville identity sector to remove');

end

