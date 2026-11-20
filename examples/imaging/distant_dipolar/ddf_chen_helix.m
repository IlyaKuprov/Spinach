% Transverse-helix frequency shifts for spherical water at 7 T, after
% Chen, Branca, and Warren, doi:10.1063/1.4721637, Figs. 1(b) and 4.
% This ideal bulk model replaces the final unwinding kick with its
% equivalent phase-sensitive receiver, omitting finite gradient lobes,
% the RF refocusing pulse, diffusion, relaxation, and radiation damping.
% The 175-degree reversal is reported in the text, not plotted in Fig. 4.
%
% The output contains evolution times (s), pulse angles (degrees), and
% dipolar-on/off complex echo ratios (time x tip x direction, axial then
% magic angle). omega contains their phase slopes (rad/s, tip x direction);
% offset_omega is the phase slope of a separate +1 Hz offset control.

function result=ddf_chen_helix()

% Spherical water and its reported dipolar time
sys.magnet=7;
sys.isotopes={'1H'};
sys.parallel={'processes',1};
sys.output='hush';
sys.disable={'hygiene'};
tau_d=0.130;
gamma=spin('1H');
magn=1/(gamma*1.25663706127e-6*tau_d);

% Spinach housekeeping
bas.formalism='sphten-liouv';
bas.approximation={'none'};
spin_system=basis(create(sys,struct()),bas);

% Encode one turn of the helix in an ideal periodic cell
p.npts=16;
p.dims=1e-3;
p.deriv={'fourier'};
p.diff=0;
p.spins={'1H'};
p.offset=0;
p.demag_factor=1/3;
p.rho_ref=state(spin_system,'Lz','1H');
p.rho0_st={p.rho_ref};
p.rho0_ph={ones(p.npts,1)};
p.rlx_ph={};
p.rlx_op={};
p.gradient_areas=[2*pi/(gamma*p.dims) 0];

% Undo the encoding in the receiver instead of applying a final kick
phase=2*pi*((0:p.npts-1)'+0.5-p.npts/2)/p.npts;
p.coil_ph={exp(1i*phase)/p.npts};
p.coil_st={coil_state(spin_system,'L+','1H')};
p.sweep=40;
p.npoints=161;
p.ddf_steps=8;
result.time=(0:p.npoints-1)'/p.sweep;

% Compare the two reported tips and a tenfold smaller tip
result.tips=[5 175 0.5];
directions=[0           0 1;...
            sqrt(2/3)   0 sqrt(1/3)];
result.echo=complex(zeros(p.npoints,numel(result.tips),size(directions,1)));
result.omega=zeros(numel(result.tips),size(directions,1));
for n=1:numel(result.tips)

    % Obtain a dipolar-off reference dividing out the preparation phase
    p.flip_angles=[deg2rad(result.tips(n)) 0];
    p.direction=directions(1,:);
    p.magnetisation=0;
    reference=imaging(spin_system,@crazed_ddf,p);

    % Compare phase evolution with the gradient on each axis
    p.magnetisation=magn;
    for k=1:size(directions,1)

        p.direction=directions(k,:);
        signal=imaging(spin_system,@crazed_ddf,p);
        result.echo(:,n,k)=signal./reference;
        coeff=polyfit(result.time,unwrap(angle(result.echo(:,n,k))),1);
        result.omega(n,k)=coeff(1);

    end

end

% Calibrate the signed phase convention against a positive Zeeman offset
p.magnetisation=0;
p.offset=1;
signal=imaging(spin_system,@crazed_ddf,p);
coeff=polyfit(result.time,unwrap(angle(signal./reference)),1);
result.offset_omega=coeff(1);

% Report signed frequency shifts with the receiver convention above
disp(table(result.tips',result.omega(:,1),result.omega(:,2),...
           'VariableNames',{'tip_deg','axial_rad_s','magic_rad_s'}));

end


