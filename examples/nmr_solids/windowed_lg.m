% Windowed LG4 for a three-proton model, using the Grenoble wlg4
% pulse phases and alternating offsets. Detection is at the end of
% each free window. The magnitude FWHM measures the dominant feature
% in the simulated spectrum, not an isolated single-spin linewidth.
%
% talos@spindynamics.org

function windowed_lg()

% Three-proton interaction model
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'1H','1H','1H'};
sys.parallel={'processes',2};
inter.zeeman.scalar={3.2 3.7 8.3};
inter.coupling.matrix=cell(3);
dipoles=[-21339.63 -95.06 56.44 97.02; -8360.57 -13.94 48.21 47.32;...
          -8939.13 2.45 52.36 -24.46];
pairs=[1 2; 1 3; 2 3];
for n=1:3
    Q=euler2dcm(deg2rad(-fliplr(dipoles(n,2:4))));
    inter.coupling.matrix{pairs(n,1),pairs(n,2)}=dipoles(n,1)*Q*diag([-1 -1 2])*Q';
end
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Rotor and frequency-switched LG settings
parameters.spins={'1H'};
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='leb_3ang_rank_5';
parameters.max_rank=4;
parameters.npoints=768;
parameters.rho0=state(spin_system,'Lx','1H');
parameters.coil=coil_state(spin_system,'L+','1H');
parameters.Ix=operator(spin_system,'Lx','1H');
parameters.Iy=operator(spin_system,'Ly','1H');
parameters.Iz=operator(spin_system,'Lz','1H');
rf_values=[80 100 120]*1e3;
windows=[3.8 4.6]*1e-6;
linewidth=zeros(numel(windows),numel(rf_values));

% Change RF and free detection window independently
for n=1:numel(windows)
    for k=1:numel(rf_values)

        % Compute the signed complex trace at the actual experimental dwell
        parameters.rf=rf_values(k);
        parameters.window=windows(n);
        dwell=2/(sqrt(1.5)*parameters.rf)+parameters.window;
        fid=singlerot(spin_system,@lg_trace,parameters,'nmr');
        assert(all(isfinite(fid))&&norm(fid)>0);

        % Measure the dominant magnitude peak between the half-height crossings
        nfft=8192;
        spectrum=abs(fftshift(fft(fid,nfft)));
        [peak,centre]=max(spectrum);
        left=find(spectrum(1:centre)<peak/2,1,'last');
        right=centre-1+find(spectrum(centre:end)<peak/2,1,'first');
        assert(~isempty(left)&&~isempty(right));
        left=left+(peak/2-spectrum(left))/(spectrum(left+1)-spectrum(left));
        right=right-(peak/2-spectrum(right))/(spectrum(right-1)-spectrum(right));
        linewidth(n,k)=(right-left)/(nfft*dwell);

    end
end

% Report linewidth at the unscaled acquisition frequencies
kfigure(); plot(rf_values/1e3,linewidth.','o-');
kxlabel('Transverse RF, kHz');
kylabel('Dominant peak magnitude FWHM, Hz'); kgrid;
klegend({'3.8 us acquisition window','4.6 us acquisition window'});

end

% Sample two opposite-offset LG pairs using midpoint rotor quadrature
function fid=lg_trace(~,parameters,H,~,~)

% Prepare the harmonic interpolation of the native rotor stack
shape=size(H{1});
samples=cellfun(@(x) full(x(:)),H,'UniformOutput',false);
coeffs=fftshift(fft([samples{:}],[],2),2)/numel(H);
rank=(numel(H)-1)/2;
harmonics=(-rank:rank).';
phase=deg2rad([55 235 125 305]);
offset=parameters.rf*[1 -1 1 -1]/sqrt(2);
rf=2*pi*parameters.rf;
pulse_time=1/(sqrt(1.5)*parameters.rf);
times=[pulse_time pulse_time parameters.window];
dwell=sum(times);
fid=complex(zeros(1,parameters.npoints));
rho=parameters.rho0;

% Propagate the finite pulses and the free window before each sample
for n=1:numel(fid)

    % Alternate the two LG pairs without rephasing the receiver
    phases=2*mod(n-1,2)+(1:2);
    drives={rf*(cos(phase(phases(1)))*parameters.Ix+...
                sin(phase(phases(1)))*parameters.Iy)+...
            2*pi*offset(phases(1))*parameters.Iz,...
            rf*(cos(phase(phases(2)))*parameters.Ix+...
                sin(phase(phases(2)))*parameters.Iy)+...
            2*pi*offset(phases(2))*parameters.Iz,0};
    elapsed=(n-1)*dwell;
    for s=1:3

        % Refine each event into midpoint slices of at most two microseconds
        slices=ceil(times(s)/2e-6);
        dt=times(s)/slices;
        for k=1:slices
            tick=elapsed+(k-0.5)*dt;
            drift=reshape(coeffs*exp(-2i*pi*parameters.rate*tick*harmonics),shape);
            U=expm(-1i*(drift+drives{s})*dt);
            rho=U*rho*U';
        end
        elapsed=elapsed+times(s);

    end
    fid(n)=trace(parameters.coil'*rho)/shape(1);

end

end


