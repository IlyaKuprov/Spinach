% Stack plotting utility for 2D NMR spectra. Syntax:
%
%   stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)
%
% Parameters:
%
%     spectrum - a real matrix containing the 2D NMR spectrum
%
%     parameters.sweep       -  one or two sweep widths, Hz
%
%     parameters.spins       -  cell array with one ot two character
%                               strings specifying the working spins
%
%     parameters.offset      -  one or two transmitter offsets, Hz
%
%     parameters.axis_units  -  axis units ('ppm','Hz','Gauss')
%
%     stack_dim              -  stacking dimension, 1 or 2
%
%     alpha_fun              -  optional function handle that takes
%                               a spectral slice and returns the al-
%                               pha value that regulates stack line
%                               opacity
%
% Outputs:
%
%     this function updates the current figure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=stack_2d.m>

function stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)

% Set common defaults
parameters=defaults(spin_system,parameters);

% Default measure of stack line opacity
if ~exist('alpha_fun','var'), alpha_fun=@(x)sqrt(norm(x,2)); end

% Check consistency
grumble(spectrum,parameters,stack_dim);

% If a complex spectrum is received, plot both components
if nnz(imag(spectrum))>0
    stack_2d(spin_system,real(spectrum),parameters,stack_dim); hold on;
    stack_2d(spin_system,imag(spectrum),parameters,stack_dim);
end

% Inform the user
report(spin_system,'plotting...');

% Accommodate homonuclear 2D sequences
if isscalar(parameters.spins)
    parameters.spins= [parameters.spins  parameters.spins];
end
if isscalar(parameters.offset)
    parameters.offset=[parameters.offset parameters.offset];
end
if isscalar(parameters.sweep)
    parameters.sweep= [parameters.sweep  parameters.sweep];
end

% Build axes and apply offsets
axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,2));
axis_f2=ft_axis(parameters.offset(2),parameters.sweep(2),size(spectrum,1));

% Convert the units 
switch parameters.axis_units
    case 'ppm'
        axis_f1=1000000*(2*pi)*axis_f1/(spin(parameters.spins{1})*spin_system.inter.magnet);
        axis_f2=1000000*(2*pi)*axis_f2/(spin(parameters.spins{2})*spin_system.inter.magnet);
        axis_f1_label=['F1: ' parameters.spins{1} ' chemical shift / ppm'];
        axis_f2_label=['F2: ' parameters.spins{2} ' chemical shift / ppm'];
    case 'Gauss'
        axis_f1=10000*(spin_system.inter.magnet-2*pi*axis_f1/spin('E'));
        axis_f2=10000*(spin_system.inter.magnet-2*pi*axis_f2/spin('E'));
        axis_f1_label='F1: magnetic induction / Gauss';
        axis_f2_label='F2: magnetic induction / Gauss';
    case 'Hz'
        axis_f1=1*axis_f1+0; axis_f2=1*axis_f2+0;
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / Hz'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / Hz'];
    case 'kHz'
        axis_f1=0.001*axis_f1+0; axis_f2=0.001*axis_f2+0;
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / kHz'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / kHz'];
    case 'MHz'
        axis_f1=0.000001*axis_f1+0; axis_f2=0.000001*axis_f2+0;
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / MHz'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / MHz'];
    case 'points'
        axis_f1=1:size(spectrum,2); axis_f2=1:size(spectrum,1);
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / points'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / points'];
    otherwise
        error('unknown axis units.');
end

% Get Cartesian arrays of points in 3D
[X,Y]=meshgrid(axis_f2,axis_f1); Z=transpose(spectrum);

% Line direction
switch stack_dim

    case 1

        % Build patch lines
        patch_x=cell([size(spectrum,2) 1]); patch_y=cell([size(spectrum,2) 1]);
        patch_z=cell([size(spectrum,2) 1]); alpha=zeros([size(spectrum,2) 1]);
        for n=1:size(Z,2)

            % Close up patch lines
            patch_x{n}=[X(:,n); NaN];
            patch_y{n}=[Y(:,n); NaN];
            patch_z{n}=[Z(:,n); NaN];

            % Get the transparency
            alpha(n)=alpha_fun(Z(:,n));

        end

    case 2

        % Build patch lines
        patch_x=cell([size(spectrum,1) 1]); patch_y=cell([size(spectrum,1) 1]);
        patch_z=cell([size(spectrum,1) 1]); alpha=zeros([size(spectrum,1) 1]);
        for n=1:size(Z,1)

            % Close up patch lines
            patch_x{n}=[X(n,:)'; NaN];
            patch_y{n}=[Y(n,:)'; NaN];
            patch_z{n}=[Z(n,:)'; NaN];

            % Get the transparency
            alpha(n)=alpha_fun(Z(n,:));

        end

end

% Normalise transparency
alpha=alpha/max(alpha); 
alpha(alpha<0.01)=0.01;

% Draw patch lines
for n=1:numel(patch_z)
    patch('XData',patch_x{n},'YData',patch_y{n},...
          'ZData',patch_z{n},'EdgeAlpha',alpha(n)); hold on;
end

% Miscellaneous cosmetics
xlim tight; ylim tight; box on; kgrid;
set(gca,'Projection','perspective'); 
camorbit(+45,-60);

% Invert the axes
set(gca,'XDir','reverse','YDir','reverse');

% Label the axes
kxlabel(axis_f2_label); kylabel(axis_f1_label);

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'axis_units')
    report(spin_system,'parameters.axis_units field not set, assuming ppm.');
    parameters.axis_units='ppm';
end
end

% Consistency enforcement
function grumble(spectrum,parameters,stack_dim)
if (~isnumeric(spectrum))||(~ismatrix(spectrum))
    error('spectrum must be a matrix.');
end
if (~isfield(parameters,'offset'))
    error('offsets should be specified in parameters.offset variable.');
end
if (numel(parameters.offset)~=1)&&(numel(parameters.offset)~=2)
    error('parameters.offset array should have one or two elements.');
end
if ~isfield(parameters,'sweep')
    error('sweep widths should be specified in parameters.sweep variable.');
end
if (numel(parameters.sweep)~=1)&&(numel(parameters.sweep)~=2)
    error('parameters.sweep array should have one or two elements.');
end
if ~isfield(parameters,'axis_units')
    error('axis units must be specified in parameters.axis_units variable.');
end
if ~ischar(parameters.axis_units)
    error('parameters.axis_units must be a character string.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if ~iscell(parameters.spins)
    error('parameters.spins should be a cell array of character strings.');
end
if (numel(parameters.spins)~=1)&&(numel(parameters.spins)~=2)
    error('parameters.spins cell array should have one or two elements.');
end
for n=1:numel(parameters.spins)
    if ~ischar(parameters.spins{n})
        error('elements of parameters.spins must be character strings.');
    end
end
if (~isnumeric(stack_dim))||(~isscalar(stack_dim))||(~ismember(stack_dim,[1 2]))
    error('stack_dim must be either 1 or 2.');
end
end

% Nobody is more inferior than those who
% insist on being equal.
%
% Friedrich Nietzsche

