% Estimates the minimum number of spatial grid points necessary to
% have a valid treatment of gradient driven experiments with expli-
% cit digitization of spatial dimensions. Syntax:
%
%            n=ngridpts(grad_amps,grad_durs,isotope,...
%                       max_coh_order,sample_size);
%
% Parameters:
%
%   grad_amps     - a row vector of all gradient amplitudes
%                   in the sequence, T/m
%
%   grad_durs     - a row vector of all gradient durations 
%                   in the sequence, s
% 
%   isotope       - the highest magnetogyric ratio isotope in
%                   the spin system, e.g. '1H'
%
%   max_coh_order - maximum order of coherence (either positive
%                   or negative) expected during the experiment
%                   being simulated
%
%   sample_size   - spatial extent of the sample, m
%
% Outputs:
%
%   n - the minimum recommended number of discretisation points
%
% Note: the function returns the minimum number of points, it may
%       in practice be necessary to have several times the number,
%       depending on your accuracy requirements.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ngridpts.m>

function n=ngridpts(grad_amps,grad_durs,isotope,max_coh_order,sample_size)

% Check consistency
grumble(grad_amps,grad_durs,isotope,max_coh_order,sample_size);

% Compute the worst-case total effective gradient
worst_case_grad=sum(abs(grad_amps.*grad_durs));

% Compute the spatial frequency of the worst-case spiral
worst_case_freq=abs(max_coh_order*spin(isotope)*worst_case_grad);

% Determine the minimum number of points required
n=ceil(worst_case_freq*sample_size/pi);

end

% Consistency enforcement
function grumble(grad_amps,grad_durs,isotope,max_coh_order,sample_size)
if (~isnumeric(grad_amps))||(~isreal(grad_amps))||(~isrow(grad_amps))
    error('grad_amps must be a row vector of real numbers.');
end
if (~isnumeric(grad_durs))||(~isreal(grad_durs))||...
   (~isrow(grad_durs))||any(grad_durs<=0)
    error('grad_durs must be a row vector of positive real numbers.');
end
if numel(grad_amps)~=numel(grad_durs)
    error('grad_amps and grad_durs must have the same number of elements.');
end
if ~ischar(isotope)
    error('isotope must be a character string.');
end
if (~isnumeric(max_coh_order))||(~isreal(max_coh_order))||...
   (~isscalar(max_coh_order))||(mod(max_coh_order,1)~=0)
    error('max_coh_order must be a real integer.');
end
if (~isnumeric(sample_size))||(~isreal(sample_size))||...
   (~isscalar(sample_size))||(sample_size<=0)
    error('sample_size must be a positive real number.');
end
end

% "At the start of every disaster movie there's a 
%  scientist being ignored."
%
% A poster brandished by an angry-looking
% scientist in Oxford's Cornmarket Street

