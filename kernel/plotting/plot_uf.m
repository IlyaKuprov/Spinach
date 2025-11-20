% Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax:
%
%             plot_uf(spin_system,spectrum_uf,parameters)
%
% Parameters:
%
%   spectrum_uf               - a real matrix containing the 2D UF NMR 
%                               spectrum
%
%   parameters.spins          - cell array with one ot two character
%                               strings specifying the working spins
%
%   parameters.dims           - sample dimension, m
%
%   parameters.deltat         - time step for the acquisition gradient, s
%
%   parameters.npoints        - number of points in the acquisition
%                               gradient
%
%   parameters.Ga             - amplitude of the acquisition gradient, T/m
%
%   parameters.offset         - two transmitter offsets for the 
%                               conventional and uf dimension, Hz                             
%
%   parameters.axis_units     - axis units ('ppm' or 'Hz')
%
%   parameters.offset_uf_cov  - offset between chemical shifts of a MQ
%                               along the F1 dimension of a conventional
%                               and an UF spectra ('ppm' or 'Hz').
%                                
% Output:
%
%   a figure with correct axis ticks axes in the UF and conventional
%   dimension
%
% mariagrazia.concilio@sjtu.edu.cn
%
% <https://spindynamics.org/wiki/index.php?title=plot_uf.m>

function plot_uf(spin_system,spectrum_uf,parameters)

% Check consistency
grumble(spectrum_uf,parameters)

% Get Ta (duration of one single acquisition gradient)
Ta=parameters.deltat*parameters.npoints; % s

% Get sweep width conventional dimension
sweep_conv=1/(2*Ta); % Hz

% Get resolution conventional dimension
res_conv=1/(2*Ta*parameters.nloops); % Hz

% Get axis in the conventional dimension
axis_f2=(-sweep_conv/2:res_conv:sweep_conv/2-1)+parameters.offset(2); 

% Get the magnetogyric ratio        
gamma=spin(parameters.spins{1}); % rad/s T

% Get kmax (maximal k-value)
k_max=gamma*parameters.Ga*Ta/(2*pi);  % m^-1

% Get constant c 
% (in according to Prog.Nucl.Magn.Reson.Spectrosc. 57,2010,241)
t_max=2*parameters.Te; 
constant_c=(-2*t_max)/parameters.dims; 

% Get resolution uf dimension
res_uf=abs(1/(constant_c*parameters.dims)); % Hz

% Get sweep width uf dimension  
sweep_uf=abs(k_max/constant_c);  % Hz

% Calculate the number of points of the Ga
ga_points_number=parameters.dims*k_max;
ga_points_number=round(ga_points_number);

% Get size of the uf dimension
uf_dim_size=size(spectrum_uf,1);

% A safety checking
if (uf_dim_size) ~= (ga_points_number)
    error('the size of the uf dimension must be equal to the product of the maximal k-value and the sample dimension.');
end

% Get axis in the uf dimension 
axis_f1=(-sweep_uf/2:res_uf:sweep_uf/2)+parameters.offset(1); 

% Get the units 
switch parameters.axis_units
    
    case 'ppm'     
          
        axis_f1=1000000*(2*pi)*axis_f1/(spin(parameters.spins{1})*spin_system.inter.magnet);
        axis_f2=1000000*(2*pi)*axis_f2/(spin(parameters.spins{2})*spin_system.inter.magnet);
        
        % Add offset between uf and conventional F1 axes 
        axis_f1=axis_f1+parameters.offset_uf_cov;
  
    case 'Hz'   
        
        axis_f1=1*axis_f1+0; axis_f2=1*axis_f2+0; 
               
        % Add offset between uf and conventional F1 axes 
        axis_f1=axis_f1+parameters.offset_uf_cov;

    otherwise
        
        % Complain and bomb out
        error('unknown axis units.');
        
end
  
% Plot the spectrum
contour(axis_f2,axis_f1,flipud(spectrum_uf));
box on; kgrid;

% Get labels
kxlabel('1Q / ppm'); kylabel('MQ / ppm');

% Invert the axes
set(gca,'XDir','reverse','YDir','reverse');

end

% Consistency enforcement
function grumble(spectrum_uf,parameters)
if (~isnumeric(spectrum_uf))||(~isreal(spectrum_uf))||(~ismatrix(spectrum_uf))
    error('the uf spectrum must be a real matrix.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if (~isfield(parameters,'offset'))
    error('offsets should be specified in parameters.offset variable.');
end
if (~isfield(parameters,'offset_uf_cov'))
    error('offset along the F1 dimension between two peaks in a conventional and UF spectra should be specified in parameters.offset_uf_cov variable.');
end
if (numel(parameters.offset)~=2)
    error('parameters.offset array should have two elements.');
end
if ~isfield(parameters,'deltat')  
    error('the time step of the acquisition gradient should be specified in the parameters.deltat variable.');  
end
if ~isfield(parameters,'npoints')  
    error('the number of points in the acquisition gradient should be specified in the parameters.npoints variable.');  
end
if ~isfield(parameters,'Ga')  
    error('the amplitude of the acquisition gradient should be specified in the parameters.Ga variable.');  
end
if ~isfield(parameters,'dims')  
    error('the sample dimension should be specified in the parameters.dims variable.');  
end
end

% Nobody has had a sabbatical at Southampton Chemistry 
% in ten years.
%
% Gill Reid (Head of Department), 
% to IK when he asked for a sabba-
% tical to write his book

