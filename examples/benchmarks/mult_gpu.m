% CPU and GPU matrix arithmetic benchmark. Set the argument to 
% either 'single' or 'double'. IK's workstation output:
%
%  GPU, precision: single
%   Titan V    PCIe (2017): 12534.7 GFLOPS
%   Tesla A100 PCIe (2021): 10285.8 GFLOPS
%   Tesla A800 PCIe (2024): 17801.2 GFLOPS
%
%  GPU, precision: double
%   Titan V    PCIe (2017): 6508.9  GFLOPS
%   Tesla A100 PCIe (2021): 10093.5 GFLOPS
%   Tesla A800 PCIe (2024): 14984.8 GFLOPS
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
sizes=power(2,20:2:26); N=sqrt(sizes);

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

