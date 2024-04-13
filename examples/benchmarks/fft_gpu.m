% GPU arithmetic benchmark - 3D Fourier transforms.
%
% i.kuprov@soton.ac.uk

function fft_gpu()

% Look for GPUs
if gpuDeviceCount==0
    disp('no CUDA GPUs detected'); return;
end

% FFT dimensions (reduce if card runs out of memory)
sizes=[128 192 256 384 512];

% Timing array
timings=zeros(numel(sizes),2,10);

% FFT size loop
for n=1:numel(sizes)
    
    % Statistics loop
    for k=1:10
        
        % CPU benchmark
        a=randn(sizes(n),sizes(n),sizes(n),'double');
        tic; fftn(a); timings(n,1,k)=toc;
        disp(['n=' num2str(sizes(n)) ', CPU time ' num2str(timings(n,1,k)) ' seconds']);
        
        % GPU benchmark
        b=randn(sizes(n),sizes(n),sizes(n),'gpuArray');
        tic; fftn(b); timings(n,2,k)=toc;
        disp(['n=' num2str(sizes(n)) ', GPU time ' num2str(timings(n,2,k)) ' seconds']);
        
    end
    
end

% Analysis
means=mean(timings(:,:,2:end),3);

% Plotting
figure(); hold on;
plot(sizes',means(:,1),'ro');
plot(sizes',means(:,2),'bo');
xlim([100 600]); box on; kgrid;
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('number of points in each dimension');
ylabel('3D FFT calculation time, seconds');
legend({'CPU','GPU'},'Location','northwest');

end
