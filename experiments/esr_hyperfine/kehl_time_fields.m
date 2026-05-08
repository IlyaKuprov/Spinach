% set the mw and rf field strength (omega1e, etc.)
% is directly called for the set up
%
% input parameters:
% constants: the Map containing the constants
% expt: the Map containing the experimental parameters
% set: boolean, should the field strength be calculated from the
%       pulse length
% values: if not fixed (set=true), array with field strengths
%
% output parameters:
% expt: updated Map with experimental parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function expt=kehl_time_fields(constants,expt,set,values)
    if nargin<4
        values=[];
    end


    % Check consistency
    grumble(constants,expt,set,values);
    % setting the field strength values as needed for the time domain ENDOR experiment
    if set==false
        if size(values,1)==2
            expt("oneE")=values(1)*2*pi*1e6;
            expt("oneN")=values(2)*2*pi*1e3;

            % omega1e in T for orientation preselection
            expt("pulsewidth")=expt("oneE")/(2*pi*constants("CONST1")*1e10);
        elseif values~=[]
            error('Wrong dimensions of input, should be: [omega_prep/2pi,omega_1e/2pi,omega_1n/2pi].');
        end
    else
        t=expt("t");
        expt("oneE")=2*pi/(4*t(1));
        if isKey(expt,"ang")
            expt("oneN")=(expt("ang")/180)*2*pi/(4*t(5));

        else
            expt("oneN")=2*pi/(4*t(5));

        end

        % omega1e in T for orientation preselection
        expt("pulsewidth")=expt("oneE")/(2*pi*constants("CONST1")*1e10);
    end
end

function grumble(constants,expt,set,values)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
if (~islogical(set))&&(~isnumeric(set))
    error('set must be logical, or numeric.');
end
if (~isempty(values))&&(~isa(values,'containers.Map'))&&(~isnumeric(values))&&(~islogical(values))
    error('values must be empty, numeric, logical, or a containers.Map object.');
end
end

