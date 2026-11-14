% Fourier normalization, Lorentzian, odd/even lengths and filter units.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_spectrum()

% Validate odd and even DFT meshes against a direct quadrature transform
result=new_test_result('kernel/cce_spectrum','Noise Fourier spectra',...
                       'Parseval, inverse normalization, Lorentzian, sidedness and units.');
p.spectrum=struct('convention','exp(+i*omega*t)','sidedness','two-sided',...
                  'window','none','frequency_grid',0,'time_domain','full','weights','rectangle');
maximum=0;
for n=[63 64]
    dt=0.03; t=((0:n-1)-floor(n/2))'*dt; c=exp(-t.^2)+0.2i*sin(2*t);
    omega=2*pi*((0:n-1)-floor(n/2))'/(n*dt); p.spectrum.frequency_grid=omega;
    s=cce_spectrum(t,c,p); dw=2*pi/(n*dt);
    time_norm=sum(abs(c).^2)*dt; frequency_norm=sum(abs(s.value).^2)*dw/(2*pi);
    result=test_close(result,'C7 numerical check',frequency_norm,time_norm,1e-12,1e-12,'Odd/even Parseval.');
    reconstructed=exp(-1i*t*omega.')*s.value*dw/(2*pi);
    result=test_close(result,'C7 numerical check',reconstructed,c,1e-12,1e-12,'Odd/even inverse angular-frequency normalization.');
    result=test_close(result,'C7 numerical check',s.frequency_hz,omega/(2*pi),0,0,'Hz axis conversion.');
    maximum=max(maximum,abs(frequency_norm-time_norm));
end

% Resolve the exponential correlation's analytic Lorentzian
rate=2; t=(0:0.0005:10)'; omega=linspace(-8,8,101)';
p.spectrum.frequency_grid=omega; p.spectrum.time_domain='positive'; p.spectrum.weights='trapezoid';
s=cce_spectrum(t,exp(-rate*t),p); reference=2*rate./(rate^2+omega.^2);
result=test_close(result,'C7 numerical check',s.value,reference,1.5e-6,1e-6,'Exponential Lorentzian quadrature.');
lorentz_error=max(abs(s.value-reference));
p.spectrum.frequency_grid=(0:0.5:8)'; p.spectrum.sidedness='one-sided';
one=cce_spectrum(t,exp(-rate*t),p);
expected=4*rate./(rate^2+p.spectrum.frequency_grid.^2); expected(1)=expected(1)/2;
result=test_close(result,'C7 numerical check',one.value,expected,3e-6,1e-6,'One-sided folding preserves the DC weight.');

% Check the zero-frequency filter integral in seconds squared
p.spectrum.modulation=ones(size(t)); filtered=cce_spectrum(t,exp(-rate*t),p);
result=test_close(result,'C7 numerical check',filtered.filter(1),100,1e-12,1e-12,'FID DC filter is total duration squared.');
refused=false;
try
    cce_spectrum(t,exp(-rate*t+1i*t),p);
catch
    refused=true;
end
result=test_true(result,'C7 contract check',refused,'Complex ordered noise cannot be folded to a one-sided PSD.');
fprintf('Spectrum Parseval error %.3g; Lorentzian max error %.3g\n',maximum,lorentz_error);

end


