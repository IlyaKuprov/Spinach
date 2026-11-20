% Solvent-referenced HOMOGENIZED spectroscopy in a linear static ramp:
% Vathyam, Lee, and Warren, Science 272, 92 (1996), DOI:
% 10.1126/science.272.5258.92. The ideal 90x-gradient-45x-echo(180x)
% sequence and the opposite indirect branches are described in Balla,
% University of Wuerzburg dissertation (2009), Fig. 3.1 and section 3.1.
%
%                    data=ddf_homogenized()
%
% Infinite bulk solvent, uniform 10 mm receiver along the ramp, complete
% grating periods across it. Water is a prescribed longitudinal field
% source after mixing; only dilute solute is detected. No diffusion,
% transverse feedback, solvent depletion, or finite RF/gradient pulses.
% The 8 Hz doublet has an unpulsed spin-half partner. The ramp runs in
% the yz plane at the magic angle to z; the correlation gradient is x.
% This geometry avoids a residual local field when the ramp winds faster
% than the encoding grating. Echo-centred direct acquisition, rotating
% frame, and fixed positive single-pulse amplitude units are used.
%
% Output: time axes (s), native FFT axes (Hz), complex single-pulse and
% 2D spectra for off/on/reversed ramps and two shift controls, isolated
% component absorptive linewidths (Hz), and half-pitches (m).

function data=ddf_homogenized()

% Set native indirect and direct bandwidths, both centred on solvent
data.t1=(0:1023)'/512;
data.t2=(0:255)/512;
data.f1=(-512:511)'/2;
data.f2=(-128:127)*2;

% Compare signed ramps (Hz/m), solvent offsets, and solute-only changes (Hz)
ramp=[0 2e4 -2e4 2e4 2e4];
common=[0 0 0 16 0];
solute=[120 120 120 120 136];
data.distance=[10 25 50 100 250 500 1000]*1e-6;
data.width=zeros(size(data.distance));
data.sq=complex(zeros(5,256));
data.spectrum=complex(zeros(1024,256,5));

% Set the solvent dipolar time, mixing angle, and echo buildup time
tau_d=1;
beta=pi/4;
echo=0.02;
decay=exp(-pi*data.t1);
direct=exp(-pi*data.t2);

% Build the native-grid spectra with the same receiver and no zero filling
for c=1:numel(ramp)

    % Integrate direct phase over the rectangular receiver along the ramp
    receive=direct.*sinc(0.01*ramp(c)*data.t2);
    offsets=solute(c)+[-4 4];
    sq=mean(exp(2i*pi*(common(c)+offsets(:))*data.t2),1).*receive;
    data.sq(c,:)=fftshift(fft(sq))/numel(data.t2);

    % Evaluate the solvent's wound Fourier mode including the Lorentz cavity
    k=pi/data.distance(1);
    winding=2*pi*ramp(c)*data.t1;
    field=1/3-(winding/sqrt(3)).^2./(k^2+winding.^2);
    angle=sin(beta)*field*(echo+data.t2)/tau_d;

    % Project the 45x-mixed, 180x-conjugated grating onto the uniform receiver
    fid=complex(zeros(size(angle)));
    for offset=offsets

        % Keep both zero-quantum branches with their hard-pulse amplitudes
        water=exp(2i*pi*common(c)*data.t1);
        reporter=exp(2i*pi*(common(c)+offset)*data.t1);
        beat=reporter.*conj(water);
        mixed=cos(beta/2)^2*conj(beat)-sin(beta/2)^2*beat;
        fid=fid+(mixed.*decay)*(exp(2i*pi*(common(c)+offset)*data.t2).*receive)/2;

    end

    % Select the first Bessel harmonic of the longitudinal solvent field
    fid=fid.*besselj(1,angle);
    data.spectrum(:,:,c)=fftshift(fft2(fid))/numel(fid);

end

% Measure the isolated strong branch at the echo, not a blended doublet
for n=1:numel(data.distance)

    % Compute the real absorption line on its acquired frequency grid
    k=pi/data.distance(n);
    winding=2*pi*ramp(2)*data.t1;
    field=1/3-(winding/sqrt(3)).^2./(k^2+winding.^2);
    line=decay.*besselj(1,sin(beta)*echo*field/tau_d);
    line(1)=line(1)/2;
    spectrum=real(fftshift(fft(line)));
    centre=find(data.f1==0);
    right=find(spectrum(centre:end)<spectrum(centre)/2,1)+centre-1;
    bins=[right-1 right];
    data.width(n)=2*interp1(spectrum(bins),data.f1(bins),spectrum(centre)/2);

end

% Present direct broadening, the resolved strong branch, and the distance sweep
if nargout==0

    % Compare on common intensity scales for the three ramp settings
    kfigure();
    tiledlayout(2,3);
    scale=max(abs(data.spectrum),[],'all');
    for c=1:3

        % Plot the direct control and solvent-referenced two-dimensional spectrum
        nexttile(c);
        plot(data.f2,abs(data.sq(c,:)));
        ylim([0 max(abs(data.sq),[],'all')]);
        kxlabel('Direct frequency, Hz');
        kylabel('SQ magnitude');
        ktitle(sprintf('Ramp: %g Hz / 10 mm',0.01*ramp(c)));
        nexttile(c+3);
        imagesc(data.f2,-data.f1,20*log10(abs(data.spectrum(:,:,c))/scale),[-60 0]);
        axis xy;
        ylim([106 134]);
        kxlabel('Direct frequency, Hz');
        kylabel('Solute minus solvent, Hz');

    end

    % Use one decibel scale relative to the strongest point in all spectra
    bar=colorbar();
    bar.Layout.Tile='east';
    bar.Label.String='Magnitude, dB';

    % Show finite-distance linewidths without imposing exact cancellation
    kfigure();
    semilogx(1e6*data.distance,data.width,'o-');
    kxlabel('Correlation distance, micrometres');
    kylabel('Component FWHM, Hz');

end

end


