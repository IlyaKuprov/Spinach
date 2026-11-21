% Test: kernel/respdor | R- and S-RESPDOR echoes
% Tests R- and S-RESPDOR against direct time-dependent propagation
% of a single 13C-14N crystallite. Syntax:
%
%                      result=test_respdor()
%
% Outputs:
%
%    result  - regression result with signed complex echo checks
%
% talos@spindynamics.org

function result=test_respdor()

% Register the reference, dephased, and dipole-off checks
result=new_test_result('kernel/respdor','R- and S-RESPDOR echoes',...
                       'Both schedules must agree with a direct Hilbert-space model.');

% Specify a tilted dipole and quadrupole at zero crystal Euler angles
sys.magnet=2*pi*600e6/spin('1H');
sys.isotopes={'13C','14N'};
inter.coupling.eigs={[] [1200 1200 -2400]; [] []};
inter.coupling.euler={[] [0.31 0.58 0]; [] []};
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{2,2}=eeqq2nqi(36e3,0,1,[-0.28 0.72 0]);
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Specify MAS timing and retain both pulse branches at each count
parameters.spins=sys.isotopes;
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.offset=[270 -750];
parameters.max_rank=24;
parameters.grid='single_crystal';
parameters.serial=true;
parameters.ncycles=[2 0 1];
parameters.sat_rf=40e3;

% Retain the active quadrupolar interaction in the dipole-off system
for n=1:2

    % Compare a coupled and an uncoupled carbon-nitrogen pair
    spin_system=test_spin_system(sys,inter,bas);
    for seq='RS'

        % Use each schedule's preparation, detection, and echo calibration
        parameters.sequence=seq;
        if seq=='R'

            % Longitudinal R echo
            parameters.echo_rf=100e3;
            parameters.sat_time=200e-6;
            parameters.rho0=state(spin_system,'Lz','13C');
            parameters.coil=-coil_state(spin_system,'Lz','13C');

        else

            % Transverse S echo
            parameters.echo_rf=80e3;
            parameters.sat_time=100e-6;
            parameters.rho0=state(spin_system,'Lx','13C');
            parameters.coil=coil_state(spin_system,'L+','13C');

        end

        % Run the public MAS path
        observed=singlerot(spin_system,@respdor,parameters,'nmr');
        if n==1

            % Compare all complex rows without using the context generator
            reference=direct_echo(parameters);
            result=test_close(result,[seq ' complex echoes'],observed,reference,2e-7,0,...
                              'Independent Hilbert-space integration, both RF branches.');
            result=test_true(result,[seq ' nonzero coupled contrast'],norm(observed(3,:))>1e-3,...
                             'Irradiation must perturb the coupled echo.');

        else

            % Require zero contrast without removing quadrupolar evolution
            result=test_close(result,[seq ' dipole-off null'],observed(3,:),zeros(1,3),1e-10,0,...
                              'A driven but uncoupled spin cannot change the carbon signal.');

        end

    end

    % Remove only the heteronuclear dipole for the null
    inter.coupling.eigs{1,2}=[];
    inter.coupling.euler{1,2}=[];

end

% Reject a window shorter than the observed-channel inversion
parameters.sat_time=0.25/parameters.echo_rf;
rejected=false;
try

    % Exercise the same validation path used by the example
    singlerot(spin_system,@respdor,parameters,'nmr');

catch err

    % Check the specific timing error rather than any exception
    rejected=contains(err.message,'parameters.sat_time must be between');

end
result=test_true(result,'reject short irradiation window',rejected,...
                 'Both branches require a non-negative flank duration.');

end

% Integrate the secular two-spin Hamiltonian in Cartesian coordinates
function curve=direct_echo(parameters)

% Construct elementary spin-1/2 and spin-1 matrices
Ix=kron([0 1; 1 0]/2,eye(3));
Iy=kron([0 -1i; 1i 0]/2,eye(3));
Iz=kron(diag([1 -1])/2,eye(3));
Sx=kron(eye(2),[0 1 0; 1 0 1; 0 1 0]/sqrt(2));
Sz=kron(eye(2),diag([1 0 -1]));

% Express the rotating tensor projections directly in Cartesian geometry
period=1/parameters.rate;
axis_z=[-parameters.axis(1),0,parameters.axis(3)];
axis_proj=@(t,a,b) axis_z*[sin(b)*cos(a-2*pi*t/period);...
                         sin(b)*sin(a-2*pi*t/period); cos(b)];
dip_hz=@(t) -1200*(3*axis_proj(t,0.31,0.58)^2-1);
quad_hz=@(t) 9e3*(3*axis_proj(t,-0.28,0.72)^2-1)/2;
offset=parameters.offset(1)*Iz+parameters.offset(2)*Sz;
drift=@(t) 2*pi*(offset+dip_hz(t)*Iz*Sz+quad_hz(t)*(3*Sz^2-2*eye(6)));

% Specify the chronological recoupling events as duration and Cartesian RF
if parameters.sequence=='R'

    % Rotary resonance and signed longitudinal receiver
    recoupling=[200e-6 10e3 0 0];
    rho=Iz;
    coil=-Iz;

else

    % R4 and its phase-inverted train at three supercycle phases
    phase=[90 270 90 270 270 90 270 90];
    phase=[phase phase+120 phase+240].';
    recoupling=[repmat(25e-6,24,1) 20e3*cosd(phase) 20e3*sind(phase) zeros(24,1)];
    rho=Ix;
    coil=Ix+1i*Iy;

end

% Build the full echo events for each irradiation arm
width=1/(2*parameters.echo_rf);
flank=(parameters.sat_time-width)/2;
delay=period-parameters.sat_time/2;
echo=[delay 0 0 0; flank 0 0 0; width parameters.echo_rf 0 0; flank 0 0 0; delay 0 0 0];
block=event_product(recoupling,drift,Ix,Iy,Sx);
echo_ref=event_product(echo,drift,Ix,Iy,Sx);
echo(2:4,4)=parameters.sat_rf;
echo_on=event_product(echo,drift,Ix,Iy,Sx);

% Evaluate direct Hilbert overlaps in the requested order
curve=repmat(trace(coil'*rho)/6,2,numel(parameters.ncycles));
for n=1:numel(parameters.ncycles)

    % Repeated rotor-synchronous arms use the same unitary on both sides
    count=parameters.ncycles(n);
    if count>0

        % Project each propagated density onto the same complex receiver
        U=block^count*echo_ref*block^count;
        V=block^count*echo_on*block^count;
        curve(:,n)=[trace(coil'*U*rho*U'); trace(coil'*V*rho*V')]/6;

    end

end
curve=[curve; curve(1,:)-curve(2,:)];

end

% Integrate each constant RF event while the drift varies continuously
function U=event_product(events,drift,Ix,Iy,Sx)

% Accumulate the Hilbert propagator across event boundaries
U=eye(6);
start=0;
options=odeset('RelTol',1e-11,'AbsTol',1e-13);
for n=1:size(events,1)

    % The R window may fill both rotor periods without an outer delay
    if events(n,1)==0
        continue;
    end

    % Integrate in seconds with RF amplitudes converted to rad/s
    stop=start+events(n,1);
    rf=2*pi*(events(n,2)*Ix+events(n,3)*Iy+events(n,4)*Sx);
    rhs=@(t,y) reshape(-1i*(drift(t)+rf)*reshape(y,6,6),36,1);
    [~,path]=ode113(rhs,[start stop],U(:),options);
    U=reshape(path(end,:),6,6);
    start=stop;

end

end


