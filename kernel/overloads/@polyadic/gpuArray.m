% Uploads all components of a polyadic object to the GPU. The object
% still looks like a polyadic to Matlab, but all of its constituent 
% matrices become gpuArrays. Syntax:
%
%                           p=gpuArray(p)
%
% Parameters:
%
%     p   - polyadic object
%
% Outputs:
%
%     p   - polyadic object with all cores,
%           prefixes and suffixes uploaded
%           to the current GPU
%
% Note: GPUs are not good at permuting array dimensions. If you find
%       yourself using polyadics, do check that the CPU isn't faster.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/gpuArray.m>

function p=gpuArray(p)

% Check consistency
grumble(p);

% Upload cores
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        p.cores{n}{k}=gpuArray(p.cores{n}{k}); 
    end
end

% Upload prefixes
for n=1:numel(p.prefix)
    p.prefix{n}=gpuArray(p.prefix{n});
end

% Upload suffixes
for n=1:numel(p.suffix)
    p.suffix{n}=gpuArray(p.suffix{n});
end

end

% Consistency enforcement
function grumble(p)
if ~isa(p,'polyadic')
    error('p must be polyadic.');
end
end

% Men are born ignorant, not stupid; they 
% are made stupid by education.
%
% Bertrand Russell

