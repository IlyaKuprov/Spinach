%NUCLEUS_GAMMA Returns legacy Kehl nuclear gyromagnetic ratios.
%
%   GAMMA=NUCLEUS_GAMMA(CONSTANTS,LABEL) returns gamma/2pi in Hz/T for
%   the isotope label used by the ENDOR port. The values preserve the
%   numerical baseline of the Kehl code while the spin quantum numbers are
%   obtained from Spinach isotope metadata elsewhere.


function gamma=kehl_nuc_gamma(constants,label)

    % Check consistency
    grumble(constants,label);
label=kehl_spin_label(label);

if strcmp(label,'1H')
    gamma=constants("GN_1H");
elseif strcmp(label,'2H')
    gamma=constants("GN_2D");
elseif strcmp(label,'14N')
    gamma=constants("GN_14N");
elseif strcmp(label,'17O')
    gamma=constants("GN_17O");
elseif strcmp(label,'19F')
    gamma=constants("GN_19F");
else
    error('unsupported nucleus label: %s',label);
end

end

function grumble(constants,label)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~ischar(label))&&(~isstring(label))
    error('label must be a character string, or a string scalar.');
end
end

