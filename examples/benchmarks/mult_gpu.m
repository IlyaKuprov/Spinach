% CPU and GPU matrix arithmetic benchmark. Set the argument to either 
% 'single' or 'double'. IK's Dell PowerEdge T630 output:
%
%  Precision: single
%  Titan V (2017, TCC mode): 753.0 GFLOPS (CPU), 12534.7 GFLOPS (GPU)
%  Tesla A100 (2021, PCI-E): 772.3 GFLOPS (CPU), 10285.8 GFLOPS (GPU)
%
%  Precision: double
%  Titan V (2017, TCC mode): 337.7 GFLOPS (CPU), 6508.9 GFLOPS (GPU)
%  Tesla A100 (2021, PCI-E): 356.4 GFLOPS (CPU), 10093.5 GFLOPS (GPU)
%
% jos.martin@mathworks.ac.uk
% ilya.kuprov@weizmann.ac.il

function mult_gpu(precision)

% Default to double precision
if ~exist('precision','var')
    precision='double';
end

% Look for GPUs
if gpuDeviceCount==0
    disp('no CUDA GPUs detected');
end

% Set matrix sizes
sizes=power(2,18:2:24); N=sqrt(sizes);

% Preallocate timing arrays
mmTimesCPU=nan(size(sizes));
mmTimesGPU=nan(size(sizes));

% Loop over sizes
for n=1:numel(sizes)
    
    % Run CPU benchmark
    A=rand(N(n),N(n),precision); 
    B=rand(N(n),N(n),precision);
    mmTimesCPU(n)=timeit(@()A*B);
    
    % Run GPU benchmark
    if gpuDeviceCount>0
        A=gpuArray(A); B=gpuArray(B);
        mmTimesGPU(n)=gputimeit(@()A*B);
    end
    
end

% Convert to GFlops
mmGFlopsCPU=(2*N.^3-N.^2)./mmTimesCPU/1e9;
mmGFlopsGPU=(2*N.^3-N.^2)./mmTimesGPU/1e9;

% Get the maxima
maxGFlopsCPU=max(mmGFlopsCPU);
maxGFlopsGPU=max(mmGFlopsGPU);

% Report to the user
disp(['Precision: ' precision]);
fprintf(['Matrix multiplication: ', ...
         '%1.1f GFLOPS (CPU), %1.1f GFLOPS (GPU)\n'], ...
         maxGFlopsCPU, maxGFlopsGPU)
      
end

