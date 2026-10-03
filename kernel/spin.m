% Literature-sourced spin properties for nuclei and magnetic particles.
% Syntax:
%
%          [gamma,multiplicity,data]=spin(name)
%
% Parameters:
%
%    name  - character row vector: nuclear label such as '13C',
%            isomer label such as '99Tc_m', particle label 'E',
%            'E+', 'N', 'M', 'M+', or 'anti1H'; '1H' is a proton
%
%            'G' is a ghost (gamma=0, multiplicity=1); 'E#' is
%            a high-spin electron with integer multiplicity #;
%            'C#', 'V#', and 'T#' are cavity, phonon, and transmon
%            modes with integer level count # and gamma=0
%
%            'table' returns the entire physical-data table in
%            the required third output; numeric outputs are empty
%
% Outputs:
%
%    gamma         - signed angular magnetogyric ratio, rad/(s*T)
%
%    multiplicity  - number of spin or population levels; 2*I+1
%                    for physical entries
%
%    data          - table row with spin, abundance, quadrupole
%                    moment, half-life (seconds), qualifiers,
%                    and field-labelled academic source DOIs;
%                    typed empty table for ghosts/abstract modes
%
% Unknown properties are missing, not zero. Numerical lookup requires a
% confirmed spin and signed gamma; other rows remain inspectable with
% spin('table'). Stable half_life=Inf means no observed decay, not proof
% of infinite lifetime. Metadata does not apply decay or abundance weights.
% The uncompressed MAT payload is loaded once per process; warm numerical
% calls use a dictionary and numeric vectors, without filesystem access.
% Clear spin or start a new MATLAB session to load an updated payload.
% Source conventions and units are in etc/isotopes_sources.md.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% a.biternas@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=spin.m>

function [gamma,multiplicity,data]=spin(name)

% Check consistency
grumble(name);

% Identify abstract modes without accessing the physical database
persistent isotopes row_map gammas mults usable
is_ghost=strcmp(name,'G');
is_mode=any(name(1)=='CVT')&&...
        ~isempty(regexp(name,'^[CVT]\d+$','once'));

% Initialise the physical table and numerical projections once per process
if isempty(isotopes)&&(~is_ghost&&~is_mode||nargout>2)
    source_file=fullfile(fileparts(mfilename('fullpath')), '..', 'etc', 'isotopes.mat');
    loaded=load(source_file,'isotopes'); isotopes=loaded.isotopes;
    if isotopes.Properties.UserData.schema~=1
        error('unsupported isotope table schema.');
    end
    row_map=dictionary(isotopes.isotope,(1:height(isotopes))');
    gammas=isotopes.gamma; mults=2*isotopes.spin+1;
    usable=isfinite(gammas)&isfinite(mults)&isotopes.spin_status=="confirmed";
end

% Return the full physical-data table only through the metadata output
if strcmp(name,'table')
    if nargout<3
        error('the table selector requires the third output.');
    end
    gamma=[]; multiplicity=[]; data=isotopes;
    return
end

% Preserve ghost and abstract mode numerical behaviour
if is_ghost
    gamma=0; multiplicity=1;
elseif is_mode
    multiplicity=str2double(name(2:end)); gamma=0;
    if multiplicity<3
        error('bosonic modes need at least three energy levels.');
    end
else

    % Resolve high-spin electron multiplicities and physical dictionary keys
    is_high=name(1)=='E'&&~isempty(regexp(name,'^E\d+$','once'));
    if is_high
        row_idx=row_map("E"); multiplicity=str2double(name(2:end));
        if multiplicity<2
            error('electrons need at least two energy levels.');
        end
    else
        row_idx=lookup(row_map,string(name),FallbackValue=0);
        if row_idx==0
            error('spin:unknown_isotope',[name ' - unknown isotope.']);
        end
        multiplicity=mults(row_idx);
    end

    % Refuse numerical simulation when required physical properties are missing
    if ~usable(row_idx)
        error('spin:data_unavailable',...
              [name ' - confirmed spin or signed magnetogyric ratio unavailable; inspect spin(''table'').']);
    end
    gamma=gammas(row_idx);

    % Construct physical metadata only when requested
    if nargout>2
        if is_high
            data=isotopes([],:);
        else
            data=isotopes(row_idx,:);
        end
    end
    return
end

% Return a typed empty metadata table for abstract specifications
if nargout>2
    data=isotopes([],:);
end

end

% Consistency enforcement
function grumble(name)
if ~ischar(name)||isempty(name)||~isrow(name)
    error('isotope specification must be a nonempty character row vector.');
end
end

% I mean that there is no way to disarm any man except through guilt.
% Through that which he himself has accepted as guilt. If a man has 
% ever stolen a dime, you can impose on him the punishment intended
% for a bank robber and he will take it. He'll bear any form of mise-
% ry, he'll feel that he deserves no better. If there's not enough 
% guilt in the world, we must create it. If we teach a man that it's
% evil to look at spring flowers and he believes us and then does it,
% we'll be able to do whatever we please with him. He won't defend 
% himself. He won't feel he's worth it. He won't fight. But save us
% from the man who lives up to his own standards. Save us from the 
% man of clean conscience. He's the man who'll beat us.
%
% Ayn Rand, "Atlas Shrugged"


