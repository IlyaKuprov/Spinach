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


function expt=kehl_cp_fields(constants,expt,set,values)
    if nargin<4
        values=[];
    end


    % Check consistency
    grumble(constants,expt,set,values);
    % setting the field strength values as needed for the CP ENDOR experiment
    t=expt("t");
    if set==false
        if size(values,2)==5
            expt("prep")=values(1)*2*pi*1e6;
            expt("SL")=values(2)*2*pi*1e6;
            expt("CP")=values(3)*2*pi*1e3;
            expt("oneE")=values(4)*2*pi*1e6;
            expt("oneN")=values(5)*2*pi*1e3;

            % omega1e in T for orientation preselection
            expt("pulsewidth")=expt("prep")/(2*pi*constants("CONST1")*1e10);
        elseif size(values,2)==3
            expt("prep")=values(1)*2*pi*1e6;
            expt("SL")=values(2)*2*pi*1e6;
            expt("CP")=values(3)*2*pi*1e3;
            expt("oneE")=2*pi/(4*t(7));
            expt("oneN")=2*pi/(2*t(5));

            % omega1e in T for orientation preselection
            expt("pulsewidth")=expt("prep")/(2*pi*constants("CONST1")*1e10);
        elseif size(values,2)==2
            if t(1)==0
                expt("prep")=2*pi/(4*t(7));
                expt("SL")=values(1)*2*pi*1e6;
                expt("CP")=values(2)*2*pi*1e3;
                % omega1e in T for orientation preselection
                expt("pulsewidth")=expt("SL")/(2*pi*constants("CONST1")*1e10);
            else
                expt("prep")=2*pi/(4*t(1));
                expt("SL")=values(1)*2*pi*1e6;
                expt("CP")=values(2)*2*pi*1e3;
                % omega1e in T for orientation preselection
                expt("pulsewidth")=expt("prep")/(2*pi*constants("CONST1")*1e10);
            end
            expt("oneE")=2*pi/(4*t(7));
            expt("oneN")=2*pi/(2*t(5));


        elseif values~=[]
            error('Wrong dimensions of input, should be: [omega_prep/2pi,omega_1e/2pi,omega_1n/2pi].');
        end
    else
        expt("prep")=2*pi/(4*t(7));
        expt("SL")=2*pi/(4*t(7));
        expt("CP")=2*pi/(2*t(5));
        expt("oneE")=2*pi/(4*t(7));
        expt("oneN")=2*pi/(2*t(5));

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

