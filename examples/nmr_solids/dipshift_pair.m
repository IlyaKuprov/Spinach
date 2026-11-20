% DIPSHIFT echoes from a 1H–13C pair at 6 kHz MAS. Compare
% finite X,Y blocks at 100 kHz RF with their matched ideal limit
% and with the distinct X-only ideal schedule. A stronger X-only
% curve illustrates the effect of changing the dipolar coupling.
%
% Calculation time: about four minutes on eight CPU workers.
%
% talos@spindynamics.org

function dipshift_pair()

% Pair specification and dipolar principal axis
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'1H','13C'};
sys.parallel={'processes',8};
S=euler2dcm(deg2rad([95.07 34.54 -90.63]));
dip_tensor=S*diag([-1 -1 2])*S';

% Basis and rotor parameters
bas.formalism='sphten-liouv';
bas.approximation={'none'};
parameters.spins=sys.isotopes;
parameters.rate=6000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=12;
parameters.grid='leb_2ang_rank_11';
parameters.time_shift=linspace(0,0.5/parameters.rate,17);

% Specify the three matched-coupling schedules and the stronger X-only case
coupling=[-10000 -10000 -10000 -23000];
width=[5e-6 0 0 0];
xy_block=[true true false false];
curves=complex(zeros(numel(coupling),numel(parameters.time_shift)));
for n=1:numel(coupling)

    % Update the system only when the physical coupling changes
    if (n==1)||(coupling(n)~=coupling(n-1))
        inter.coupling.matrix=cell(2,2);
        inter.coupling.matrix{1,2}=coupling(n)*dip_tensor;
        spin_system=create(sys,inter);
        spin_system=basis(spin_system,bas);
        parameters.rho0=state(spin_system,'Lx','13C');
        parameters.coil=coil_state(spin_system,'L+','13C');
    end

    % Propagate and scale by the prepared carbon overlap, not by the echo
    parameters.pulse_width=width(n);
    parameters.xy_block=xy_block(n);
    signal=singlerot(spin_system,@dipshift,parameters,'nmr');
    curves(n,:)=signal/(parameters.coil'*parameters.rho0);

end

% Plot the signed echo against fractional rotor displacement
kfigure();
plot(parameters.rate*parameters.time_shift,real(curves));
kxlabel('Pulse displacement, rotor periods');
kylabel('Signed echo / initial overlap');
klegend({'Finite X,Y, -10 kHz','Ideal X,Y, -10 kHz',...
         'Ideal X-only, -10 kHz','Ideal X-only, -23 kHz'});
kgrid; xlim tight;

end


