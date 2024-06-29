% Returns a cell array of drift Liouvillians suitable for the
% control.drifts variable in ensemble control optimisations.
% Syntax:
%
%  [drifts,spc_dim]=drifts(spin_system,context,parameters)
%
% Parameters:
%
%     context     - a function handle to Spinach context
%                   responsible for handling the ensemble
%                   (@powder, @singlerot, etc.)
%
%     parameters  - parameters required by the context
%
%     assumptions - assumptions required by the context
%
% Outputs:
%
%     drifts      - a cell array of Liouvillians format-
%                   ted as {{La},{Lb},...}, one per en-
%                   semble member
%
%     spc_dim     - dimension of the classical dynamics
%                   subspace (e.g. rotor grid in MAS)
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=drifts.m>

function [drifts,spc_dim]=drifts(spin_system,context,...
                                 parameters,assumptions)
                             
% Check consistency
grumble(context,assumptions);

% Disable ensemble summation
parameters.sum_up=0;
parameters.serial=1;

% Capture evolution generators
gen_cap=@(varargin)varargin(3:end);
                             
% Call the user-specified context
systems=context(spin_system,gen_cap,parameters,assumptions);

% Pull drift Liouvillians
drifts=cell(1,numel(systems));
for n=1:numel(systems)
    
    % Get Liouvillian components
    H=systems{n}{1}; R=systems{n}{2}; K=systems{n}{3};
    
    % Assign drift Liouvillians
    drifts{n}={H+1i*R+1i*K};
    
    % Get hydrodynamics if present
    if numel(systems{n})==5
        drifts{n}{1}=drifts{n}{1}+1i*systems{n}{5};
    end
    
end

% Return classical subspace dimension
spc_dim=size(H,1)/size(spin_system.bas.basis,1);

end

% Consistency enforcement
function grumble(context,assumptions)
if ~isa(context,'function_handle')
    error('context must be a function handle.');
end
if ~ischar(assumptions)
    error('assumptions must be a character string.');
end
end

% And so first she was crying for a long time, and
% then she became evil.
%
% Mikhail Bulgakov, "Master and Margarita"

