% Resolves the complex line spectrum of a tick-sampled unitary period.
% The convention is signal(t)=sum(amplitude.*exp(-2i*pi*frequency*t)).
% Frequencies are in Hz, folded into [-N*rate/2,N*rate/2), where N
% is the number of ticks per period and rate is the positive rotor
% frequency. Sideband indices refer to the folded secular frequencies
% in [-rate/2,rate/2); the total frequency may wrap at the tick Nyquist.
%
% Syntax:
%
%   [frequency,amplitude,sideband]=period_lines(prefix,rate,rho,coil)
%
% Parameters:
%
%   prefix - row cell of equally spaced unitary prefixes, including
%            identity and the complete period
%   rate   - positive rotor rate in Hz
%   rho    - initial density matrix
%   coil   - detection matrix, not necessarily Hermitian
%
% Outputs: equally sized columns; coincidences are summed before
% returning, with no amplitude threshold.
%
% talos@spindynamics.org

function [frequency,amplitude,sideband]=period_lines(prefix,rate,rho,coil)

% Check consistency
grumble(prefix,rate,rho,coil);

% Transform density and time-dependent detection into the period basis
[V,T]=schur(full(prefix{end}),'complex');
phase=angle(diag(T));
rho=V'*rho*V;
count=numel(prefix)-1;
tick=(0:count-1)/(count*rate);
dim=size(rho,1);
detection=complex(zeros(dim^2,count));
for n=1:count
    U=prefix{n}*V;
    detector=U'*coil'*U;
    detection(:,n)=reshape(detector.',[],1);
end

% Remove the folded secular phase before Fourier transforming micromotion
base=reshape(-(phase-phase.')*rate/(2*pi),[],1);
base=mod(base+rate/2,rate)-rate/2;
local=detection.*rho(:).*exp(2i*pi*base*tick);
amplitude=fftshift(ifft(local,[],2),2);
sideband=repmat(-floor(count/2):ceil(count/2)-1,dim^2,1);
bandwidth=count*rate;
frequency=mod(base+rate*sideband+bandwidth/2,bandwidth)-bandwidth/2;

% Combine coincidences at the scale of angular roundoff
[frequency,order]=sort(frequency(:));
tolerance=32*eps(bandwidth);
group=cumsum([true;diff(frequency)>tolerance]);
first=[true;diff(group)>0];
amplitude=accumarray(group,amplitude(order));
sideband=sideband(order);
sideband=sideband(first);
frequency=frequency(first);

end

% Consistency enforcement
function grumble(prefix,rate,rho,coil)
validateattributes(rate,{'numeric'},{'scalar','real','finite','positive'});
if ~iscell(prefix)||~isrow(prefix)||numel(prefix)<2
    error('Supply the identity and all propagated prefixes in a cell row.');
end
if ~isnumeric(rho)||~isnumeric(coil)||~ismatrix(rho)||...
   size(rho,1)~=size(rho,2)||~isequal(size(rho),size(coil))
    error('The initial state and detection matrix must have the same square dimensions.');
end
dim=size(rho,1);
if ~all(cellfun(@(x) isnumeric(x)&&isequal(size(x),[dim dim])&&...
    norm(x'*x-eye(dim),'fro')<=32*numel(prefix)*dim*eps,prefix))
    error('The prefixes must be unitary matrices of the density matrix dimension.');
end
if norm(prefix{1}-eye(dim),'fro')>32*dim*eps
    error('The first prefix must be the identity.');
end
end


