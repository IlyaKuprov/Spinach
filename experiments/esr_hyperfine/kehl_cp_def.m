% set the CP dimension
% is directly called for the set up
%
% input parameters:
% constants: the Map containing the constants
% expIN: the Map containing the experimental parameters
% start_CP: start frequency for CP pulse
% Npts_CP: number of points for CP axis
% range_CP: sweep range of CP axis
%
% output parameters:
% expt: updated Map with experimental parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function expt=kehl_cp_def(expIN,start_CP,Npts_CP,range_CP)

    % Check consistency
    grumble(expIN,start_CP,Npts_CP,range_CP);
    expt=expIN;
    expt("start_CP")=start_CP*1e6;
    expt("Npts_CP")=Npts_CP;
    expt("range_CP")=range_CP*1e6;
end

function grumble(expIN,start_CP,Npts_CP,range_CP)
if ~isa(expIN,'containers.Map')
    error('expIN must be a containers.Map object.');
end
if ~isnumeric(start_CP)
    error('start_CP must be numeric.');
end
if ~isnumeric(Npts_CP)
    error('Npts_CP must be numeric.');
end
if ~isnumeric(range_CP)
    error('range_CP must be numeric.');
end
end

