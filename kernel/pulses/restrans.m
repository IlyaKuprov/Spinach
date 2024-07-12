% RLC circuit response calculation - converts a waveform from the
% ideal shape emitted by the instrument into the shape that comes
% out of the RLC circuit of the probe. Syntax:
%
%          [X,Y,dt]=restrans(X_user,Y_user,dt_user,...
%                            omega,Q,model,up_factor)
%
% Parameters:
%
%       X_user    - in-phase part of the rotating frame
%                   pulse waveform, a column vector of 
%                   real numbers
%
%       Y_user    - out-of-phase part of the rotating 
%                   frame pulse waveform, a column vec-
%                   tor of real numbers
%
%       dt_user   - time slice duration, seconds
%
%       omega     - RLC circuit resonance frequency in
%                   radians per second, a real number
%
%       Q         - RLC circuit quality factor, a real
%                   positive number
%
%       model     - input signal model, use 'pwc' for
%                   piecewise-constant, and 'pwl' for
%                   piecewise-linear input; time shift
%                   compensation for piecewise-linear
%                   is requested by 'pwl_tsc'
%
%       up_factor - the output waveform will have more
%                   discretisation points than the in-
%                   put waveform by this factor, about
%                   100 is a safe guess
%
% Outputs:
%
%       X        - in-phase part of the rotating frame
%                  pulse waveform distorted by the RLC
%                  response, a column vector of real 
%                  numbers
%
%       Y        - out-of-phase part of the rotating
%                  frame pulse waveform distorted by
%                  the RLC response, a column vector
%                  of real numbers
%
%       dt       - slice duration in the distorted wave-
%                  form, seconds
%
% u.rasulov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=restrans.m>

function [X,Y,dt]=restrans(X_user,Y_user,dt_user,omega,Q,model,up_factor)

% Check consistency
grumble(X_user,Y_user,dt_user,omega,Q)

% 16 x Nyquist oversampling
dt=pi/(16*omega);

% Signal model
switch model

    % Piecewise-constant
    case 'pwc'

        % Assume that the user has slice midpoints
        nslices=numel(X_user); tmax=nslices*dt_user;
        user_time_grid=dt_user*(cumsum(ones(nslices,1))-1/2);
        
        % Time grid required by the RLC circuit
        circuit_time_grid=linspace(0,tmax,tmax/dt+1)';

        % Interpolate user waveform onto RLC circuit time grid
        X0=interp1(user_time_grid,X_user,circuit_time_grid,'nearest','extrap');
        Y0=interp1(user_time_grid,Y_user,circuit_time_grid,'nearest','extrap');
        
    % Piecewise-linear
    case {'pwl','pwl_tsc'}

        % Assume that the user has slice edges
        nslices=numel(X_user)-1; tmax=nslices*dt_user;
        user_time_grid=linspace(0,tmax,nslices+1);
        
        % Time grid required by the RLC circuit
        circuit_time_grid=linspace(0,tmax,tmax/dt+1)';

        % Interpolate user waveform onto RLC circuit time grid
        X0=interp1(user_time_grid,X_user,circuit_time_grid,'linear');
        Y0=interp1(user_time_grid,Y_user,circuit_time_grid,'linear');

    otherwise

        % Complain and bomb out
        error('unknown signal model.');

end

% Fancy colour palette
rgb_rd=[0.8500 0.3250 0.0980]; 
rgb_bd=[0.0000 0.4470 0.7410];
rgb_rb=rgb_rd+1; rgb_rb=rgb_rb/max(rgb_rb);
rgb_bb=rgb_bd+1; rgb_bb=rgb_bb/max(rgb_bb);

% Diagnostic plotting if no outputs
if nargout==0
    plot(user_time_grid,X_user,'.','Color',rgb_rb); hold on;
    plot(user_time_grid,Y_user,'.','Color',rgb_bb);
    plot(circuit_time_grid,X0,'-','Color',rgb_rb); 
    plot(circuit_time_grid,Y0,'-','Color',rgb_bb); 
end

% Generate wall clock input signal
inp_amp=sqrt(X0.^2+Y0.^2); inp_phi=atan2(Y0,X0);
inp_signal=inp_amp.*cos(omega*circuit_time_grid+inp_phi);

% Build the RLC circuit response kernel
sys=tf(1/Q,[1/(omega^2) 1/(omega*Q) 1]);

% Apply the RLC circuit response kernel
out_signal=lsim(sys,inp_signal,circuit_time_grid);

% Heterodyne out the carrier frequency
X=2*lowpass(out_signal.*sin(omega*circuit_time_grid),...
            1,64,ImpulseResponse="iir",Steepness=0.95);
Y=2*lowpass(out_signal.*cos(omega*circuit_time_grid),...
            1,64,ImpulseResponse="iir",Steepness=0.95);

% Downsample the rotating frame waveform
down_factor=numel(X)/(up_factor*numel(X_user));
down_factor=floor(down_factor); dt=dt*down_factor;
X=X(1:down_factor:end); Y=Y(1:down_factor:end);
circuit_time_grid=circuit_time_grid(1:down_factor:end);

% Do time shift compensation
if strcmp(model,'pwl_tsc')
    circuit_time_grid=circuit_time_grid-2*Q/omega;
end

% Diagnostic plotting if no outputs
if nargout==0
    plot(circuit_time_grid,X,'-','Color',rgb_rd);
    plot(circuit_time_grid,Y,'-','Color',rgb_bd);
    kxlabel('time, seconds'); xlim('tight');
    kylabel('amplitude, rad/s'); kgrid;
    legend({'X, user',  'Y, user', ...
            'X, input', 'Y, input',...
            'X, output','Y, output'},  ...
            'Location','south');
end

end

% Consistency enforcement
function grumble(X_user,Y_user,dt_user,omega,Q)
if (~isnumeric(dt_user))||(~isreal(dt_user))||...
   (~isscalar(dt_user))||(dt_user<=0)
    error('dt_user must be a real positive scalar.');
end
if (~isnumeric(X_user))||(~isreal(X_user))||(~iscolumn(X_user))||...
   (~isnumeric(Y_user))||(~isreal(Y_user))||(~iscolumn(Y_user))
    error('X_user and Y_user must be real column vectors.');
end
if (~isnumeric(omega))||(~isreal(omega))||(~isscalar(omega))
    error('omega must be a real scalar.');
end
if (~isnumeric(Q))||(~isreal(Q))||(~isscalar(Q))||(Q<=0)
    error('Q must be a positive real scalar.');
end
if dt_user<(pi/omega)
    error('dt_user breaks rotating frame approximation.');
end
end

% Criminals are a small minority in any age or country. And the harm
% they have done to mankind is infinitesimal compared to the horrors -
% the bloodshed, the wars, the persecutions, the confiscations, the 
% famines, the enslavements - perpetrated by mankind's governments.
%
% Ayn Rand

