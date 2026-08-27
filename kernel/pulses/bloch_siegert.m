% Applies Bloch-Siegert corrections to Cartesian control pulses. Takes
% control operators and their amplitude arrays and returns augmented
% arrays with the Bloch-Siegert response operator of each control chan-
% nel appended as a virtual channel whose coefficient vector is the
% square of the corresponding control amplitude vector. The output
% arguments are ready for the shaped_pulse_xy function, which then
% propagates the system under the same slice generators that the GRAPE
% engines use when Bloch-Siegert corrections are enabled. Syntax:
%
%       [ctrl_opers,...
%        ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers,...
%                                              ctrl_coefs)
%
% Parameters:
%
%    spin_system - Spinach spin system object containing the
%                  Bloch-Siegert settings added by optimcon()
%
%     ctrl_opers - a cell array of control operators, one per
%                  control channel
%
%     ctrl_coefs - a cell array of control coefficient vectors
%                  in rad/s, one vector per control channel
%
% Outputs:
%
%     ctrl_opers - an augmented cell array of control opera-
%                  tors, with the response operator of each
%                  channel appended
%
%     ctrl_coefs - an augmented cell array of coefficient vec-
%                  tors, with the squared amplitude vector of
%                  each channel appended
%
% Note: the augmented arrays are only valid for the piecewise-constant
%       propagation methods of shaped_pulse_xy; piecewise-linear methods
%       would interpolate the squared coefficients, which does not cor-
%       respond to the square of the interpolated control amplitude.
%
% aditya.dev@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bloch_siegert.m>

function [ctrl_opers,...
          ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers,...
                                                ctrl_coefs)

% Check consistency
grumble(spin_system,ctrl_opers,ctrl_coefs);

% Count the physical channels
n_ctrls=numel(ctrl_opers);

% Append the response operator channels
for n=1:n_ctrls
    ctrl_opers{n_ctrls+n}=spin_system.control.resp_ops{n};
    ctrl_coefs{n_ctrls+n}=ctrl_coefs{n}.^2;
end

end

% Consistency enforcement
function grumble(spin_system,ctrl_opers,ctrl_coefs)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if ~isfield(spin_system,'control')
    error('control data missing from spin_system, run optimcon() first.');
end
if (~isfield(spin_system.control,'bsiegert'))||...
   (~islogical(spin_system.control.bsiegert))||...
   (~isscalar(spin_system.control.bsiegert))||...
   (~spin_system.control.bsiegert)
    error('Bloch-Siegert corrections must be enabled in optimcon().');
end
if (~isfield(spin_system.control,'resp_ops'))||...
   (~iscell(spin_system.control.resp_ops))
    error('response operators missing, run optimcon() first.');
end
if (~iscell(ctrl_opers))||isempty(ctrl_opers)||...
   (numel(ctrl_opers)~=numel(spin_system.control.resp_ops))
    error('ctrl_opers must be a cell array with one element per control channel.');
end
if (~iscell(ctrl_coefs))||(numel(ctrl_coefs)~=numel(ctrl_opers))
    error('ctrl_coefs must be a cell array with one element per control channel.');
end
for n=1:numel(ctrl_opers)
    if (~isnumeric(ctrl_opers{n}))||...
       (size(ctrl_opers{n},1)~=size(ctrl_opers{n},2))
        error('all elements of ctrl_opers must be square matrices.');
    end
    if (~isnumeric(ctrl_coefs{n}))||(~isreal(ctrl_coefs{n}))||...
       (~isvector(ctrl_coefs{n}))
        error('all elements of ctrl_coefs must be real vectors.');
    end
    if numel(ctrl_coefs{n})~=numel(ctrl_coefs{1})
        error('all elements of ctrl_coefs must have the same length.');
    end
end
end

% During a state visit to India, the Shah of Iran and
% the King of Nepal were overheard talking over dinner
% about Iran's mounting tensions. The King of Nepal
% asked the Shah: "I hear the mullahs and maulawis in
% Iran speak a great deal against your regime. Is it
% true?" - "Yes. Things are... tense. But not all of
% them speak against us. Some are very quiet." - "Who
% are the quiet ones?" The Shah lifts his glass: "The
% ones in their graves."
%
% Story courtesy of AD's grandfather

