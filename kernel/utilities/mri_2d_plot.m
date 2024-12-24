% 2D MRI image plotting. Syntax: 
%
%                 mri_2d_plot(mri,parameters,method)
% 
% Parameters: 
%
%      method            - 'image' uses gradient information to
%                           determine field of view; 'phantom' uses
%                           phantom dimensions; 'k-space' assumes
%                           that a k-space representation has been
%                           supplied and uses gradient information
%                           and the real part is plotted
%
%      mri               - 2D MRI image or a phantom, or the k-space
%                          representation of a 2D MRI image
%
%      parameters.spins  - nuclei on which the sequence ran.
%
%      parameters.pe_grad_amp  - amplitude of the phase encoding
%                                gradient, T/m
%
%      parameters.pe_grad_dur  - duration of the phase encoding
%                                gradient, seconds
%
%      parameters.ro_grad_amp  - the amplitude of the readout
%                                gradient, T/m
%
%      parameters.ro_grad_dur  - the duration of the readout
%                                gradient, seconds.
%
% m.g.concilio@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mri_2d_plot.m>

function mri_2d_plot(mri,parameters,method)

% Check consistency
grumble(mri,parameters,method)

% Set the option
switch method
    
    case 'image'
        
        % Get the magnetogyric ratio
        gamma=spin(parameters.spins{1});
        
        % Get field of view along dimension 1 (phase encoding)
        max_freq=gamma*parameters.pe_grad_dur*parameters.pe_grad_amp;
        pixel_size=1*pi/max_freq; fov1=size(mri,1)*pixel_size;

        % Get field of view along dimension 2 (frequency encoding)
        max_freq=gamma*parameters.ro_grad_dur*parameters.ro_grad_amp;
        pixel_size=2*pi/max_freq; fov2=size(mri,2)*pixel_size;

        % Compute axis extents
        fov1_axis=[-fov1/2,+fov1/2]; fov2_axis=[-fov2/2,+fov2/2];

        % Do the plotting
        imagesc(fov2_axis,fov1_axis,mri); 

    case 'phantom'

        % Get axis extents
        fov1_axis=[-parameters.dims(1)/2,+parameters.dims(1)/2];
        fov2_axis=[-parameters.dims(2)/2,+parameters.dims(2)/2];

        % Do the plotting
        imagesc(fov2_axis,fov1_axis,mri); 

    case 'k-space'

        % Get the magnetogyric ratio
        gamma=spin(parameters.spins{1});
        
        % Get max spatial frequency along dimension 1 (phase encoding)
        max_freq_1=gamma*parameters.pe_grad_dur*parameters.pe_grad_amp;
        
        % Get max spatial frequency along dimension 2 (frequency encoding)
        max_freq_2=gamma*parameters.ro_grad_dur*parameters.ro_grad_amp;
        
        % Compute axis extents in Hz/m
        k1_axis=[-max_freq_1/2,+max_freq_1/2]/(2*pi); 
        k2_axis=[-max_freq_2/2,+max_freq_2/2]/(2*pi);

        % Plot the real part of the k-space signal
        imagesc(k2_axis,k1_axis,real(mri)); 

    otherwise   
        
        % Complain and bomb out  
        error('unknown plotting method.');

end

% Do axis labels
switch method

    case {'image','phantom'}

        % Fields of view in metres
        kxlabel('FOV2, m'); kylabel('FOV1, m');

    case 'k-space'

        % Spatial frequency ranges, Hz/m
        kxlabel('$k_2$, rad/m'); kylabel('$k_1$, rad/m');

end

% Apply the colour map 
colormap(contrast([0 1])); 

% Tidy up axis limits
axis equal; xlim tight; 
ylim tight; drawnow();

end

% Consistency enforcement
function grumble(mri,parameters,method)
if ~ischar(method), error('method must be a character string.'); end
if ismember(method,{'image','phantom'})
    if (~isnumeric(mri))||(~isreal(mri))||(~ismatrix(mri))
        error('mri must be a real matrix.');
    end
end
if ismember(method,{'image','k-space'})
    if ~isfield(parameters,'ro_grad_amp')
        error('readout gradient amplitude should be specified in parameters.ro_grad_amp variable.');
    end
    if ~isfield(parameters,'ro_grad_dur')
        error('readout gradient duration should be specified in parameters.ro_grad_dur variable.');
    end
    if ~isfield(parameters,'pe_grad_amp')
        error('phase encoding gradient amplitude should be specified in parameters.pe_grad_amp variable.');
    end
    if ~isfield(parameters,'pe_grad_dur') 
        error('phase encoding duration gradient should be specified in parameters.pe_grad_dur variable.');
    end
    if numel(parameters.ro_grad_amp)~=1
        error('parameters.ro_grad_amp should have exactly one element.');
    end    
    if numel(parameters.ro_grad_dur)~=1
        error('parameters.ro_grad_dur should have exactly one element.');
    end   
    if numel(parameters.pe_grad_amp)~=1 
        error('parameters.pe_grad_amp should have exactly one element.');
    end
    if numel(parameters.pe_grad_dur)~=1  
        error('parameters.pe_grad_dur should have exactly one element.');
    end
end
if strcmp(method,'phantom')
    if ~isfield(parameters,'dims')
        error('phantom dimensions must be specified in parameters.dims variable.');
    end  
end
end

% "But I am very poorly today and very stupid and hate 
% everybody and everything."
%
% From Charles Darwin's diaries

