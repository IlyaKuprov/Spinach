% GPU communications benchmark. Adapted from example code in Matlab
% documentation. Data for IK's favourite NVIDIA cards on Dell Power
% Edge T630 workstation:
%
% Tesla K40  (2013):           send 9.7,  gather 2.6, bw GPU 190,  bw CPU 64
% Titan V    (2017, TCC mode): send 10.3, gather 2.6, bw GPU 568,  bw CPU 64
% Tesla A100 (2021, PCI-E):    send 10.4, gather 2.6, bw GPU 1291, bw CPU 64
%
% i.kuprov@soton.ac.uk

function comm_gpu(n)

% All GPUs by default
if nargin==0
    for n=1:gpuDeviceCount('available')
        comm_gpu(n);
    end
    return
end

% Pick the GPU
gpu=gpuDevice(n); 
namestring=['GPU ' num2str(n) ': ' gpu.Name];
disp(namestring);

% 8 bytes per double
size_of_double=8;

% Array sizes to test
sizes=power(2,14:28);

% Preallocate answer arrays
send_times=inf(size(sizes));
gather_times=inf(size(sizes));
memory_times_gpu=inf(size(sizes));
memory_times_host=inf(size(sizes));

% Measure performance
for n=1:numel(sizes)
    
    % Generate arrays
    num_elements=sizes(n)/size_of_double;
    host_data=randi([0 9],num_elements,1);
    gpu_data=randi([0 9],num_elements,1,'gpuArray');
    
    % Time sending to GPU
    send_fcn=@()gpuArray(host_data);
    send_times(n)=gputimeit(send_fcn);
    
    % Time gathering from GPU
    gather_fcn=@()gather(gpu_data);
    gather_times(n)=gputimeit(gather_fcn);
    
    % Time read+write inside GPU
    plus_fcn_gpu=@()plus(gpu_data,1.0);
    memory_times_gpu(n)=gputimeit(plus_fcn_gpu);
    
    % Time read+write on host
    plus_fcn_host=@()plus(host_data,1.0);
    memory_times_host(n)=timeit(plus_fcn_host);
    
end

% Get performance figures
send_bandwidth=(sizes./send_times)/1e9;
max_send_bandwidth=max(send_bandwidth);
gather_bandwidth=(sizes./gather_times)/1e9;
max_gather_bandwidth=max(gather_bandwidth);
memory_bandwidth_gpu=2*(sizes./memory_times_gpu)/1e9;
max_bwgpu=max(memory_bandwidth_gpu);
memory_bandwidth_host=2*(sizes./memory_times_host)/1e9;
max_bwhost=max(memory_bandwidth_host);

% Report to console
disp(['Peak send bandwidth:       ' num2str(max_send_bandwidth) ' GB/s']);
disp(['Peak gather bandwidth:     ' num2str(max_gather_bandwidth) ' GB/s']);
disp(['Peak RW bandwidth on GPU:  ' num2str(max_bwgpu) ' GB/s']);
disp(['Peak RW bandwidth on host: ' num2str(max_bwhost) ' GB/s']);

% Do the plotting
figure('Name',namestring,'NumberTitle','off'); 
subplot(1,2,1); hold on; kgrid; box on;
plot(sizes,send_bandwidth,'bo',...
     sizes,gather_bandwidth,'ro');
legend({'To GPU','From GPU'},'Location','NorthWest');
set(gca,'XScale','log'); xlabel('Array size (bytes)');
ylabel('Transfer bandwidth (GB/s)');
subplot(1,2,2); hold on; kgrid; box on;
plot(sizes,memory_bandwidth_gpu,'bo',...
     sizes,memory_bandwidth_host,'ro');
legend('On GPU','On Host','Location','NorthWest')
set(gca,'XScale','log','YScale','log'); 
xlabel('Array size (bytes)'); ylabel('Bandwidth (GB/s)'); 
scale_figure([1.5 0.75]); drawnow();

end

