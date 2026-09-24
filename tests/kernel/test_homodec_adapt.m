% Tests acquisition irradiation through Hilbert admission. Syntax:
%
%                   result=test_homodec_adapt()
%
% Inputs: none
%
% Outputs:
%
%    result - regression checks for a noncommuting 1H-13C system
%
% A complex-phase soft pulse and transverse acquisition irradiation
% must give the same FID after Hilbert admission and native Liouville
% execution. Zero power, absent irradiation, and no-op formalisms are
% checked separately. Optional operator caching is not enabled.
%
% talos@spindynamics.org

function result=test_homodec_adapt()

% Initialise the regression record
result=new_test_result('kernel/homodec_adapt',...
                       'Acquisition irradiation admission',...
                       'Hilbert irradiation must become a commutation superoperator.');

% Build coupled nuclei with unequal offsets
sys.magnet=9.4;
sys.isotopes={'1H','13C'};
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'processes',1};
inter.zeeman.scalar={1 3};
inter.coupling.scalar={0 150;0 0};
bas.approximation='none';
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv','zeeman-wavef'};
fids=cell(2,3);

% Compare admitted and native requests with identical physical parameters
for n=1:numel(formalisms)

    % Build operators in the specified representation
    bas.formalism=formalisms{n};
    spin_system=basis(create(sys,inter),bas);
    spin_system=assume(spin_system,'nmr');
    H=hamiltonian(spin_system);
    parameters.homodec_oper=cos(0.4)*operator(spin_system,'Lx','1H')+...
                           sin(0.4)*operator(spin_system,'Ly','1H');
    parameters.homodec_pwr=100;

    % Require genuinely noncommuting drift and irradiation
    comm_norm=norm(full(H*parameters.homodec_oper-parameters.homodec_oper*H),'fro');
    result=test_true(result,'noncommuting irradiation',comm_norm>1,...
                     'offsets and scalar coupling do not commute with transverse irradiation');

    % Check conversion and unchanged native formalism inputs
    [converted,admitted]=sim2liouv(spin_system,parameters,H,[],[]);
    if n==1
        expected=hilb2liouv(parameters.homodec_oper,'comm');
        result=test_close(result,'commutation operator',admitted.homodec_oper,...
                          expected,1e-14,1e-14,'irradiation uses the standard commutator conversion');
    else
        result=test_true(result,'no-op formalism',...
                         isequaln(converted,spin_system)&&isequaln(admitted,parameters),...
                         'native Liouville and wavefunction inputs are returned unchanged');
    end
    if n>2, continue; end

    % Specify a soft pulse and non-Hermitian detection operator
    parameters.spins={'1H'};
    parameters.offset=100;
    parameters.sweep=2000;
    parameters.npoints=8;
    parameters.rho0=state(spin_system,'Lz','1H');
    parameters.coil=state(spin_system,'L+','1H');
    parameters.pulse_frq=50;
    parameters.pulse_phi=0.3;
    parameters.pulse_pwr=2*pi*1000;
    parameters.pulse_dur=0.00025;
    parameters.pulse_rnk=2;
    parameters.method='expm';
    parameters.dead_time=0.00013;
    fids{n,1}=liquid(spin_system,@sp_acquire,parameters,'nmr');

    % Compare zero irradiation with complete absence of the optional fields
    parameters.homodec_pwr=0;
    fids{n,2}=liquid(spin_system,@sp_acquire,parameters,'nmr');
    parameters=rmfield(parameters,{'homodec_oper','homodec_pwr'});
    fids{n,3}=liquid(spin_system,@sp_acquire,parameters,'nmr');
    [~,admitted]=sim2liouv(spin_system,parameters,H,[],[]);
    result=test_true(result,'absent operator',~isfield(admitted,'homodec_oper'),...
                     'admission does not invent optional irradiation');
    result=test_close(result,'zero power',fids{n,2},fids{n,3},1e-10,1e-10,...
                      'zero irradiation gives the unirradiated signal');
end

% Check physical signal and a measurable acquisition-stage irradiation effect
result=test_true(result,'nonzero signal',norm(fids{2,1})>0.1,...
                 'formalism agreement is not a comparison of vanishing signals');
result=test_true(result,'irradiation effect',norm(fids{2,1}-fids{2,3})>0.01,...
                 'the acquisition-stage field changes the measured signal');
for n=1:3
    result=test_close(result,'Hilbert/native FID',fids{1,n},fids{2,n},1e-7,1e-7,...
                      'admitted and native Liouville trajectories agree');
end

end


