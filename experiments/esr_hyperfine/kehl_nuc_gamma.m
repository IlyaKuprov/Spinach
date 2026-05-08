%NUCLEUS_GAMMA Returns legacy Kehl nuclear gyromagnetic ratios.
%
%   GAMMA=NUCLEUS_GAMMA(LABEL) returns gamma/2pi in Hz/T for the isotope
%   label used by the ENDOR port. The values preserve the numerical
%   baseline of the Kehl code while the spin quantum numbers are obtained
%   from Spinach isotope metadata elsewhere.


function gamma=kehl_nuc_gamma(label)

    % Check consistency
    grumble(label);
label=kehl_spin_label(label);

if strcmp(label,'1H')
    gamma=42.576E6;
elseif strcmp(label,'2H')
    gamma=6.536E6;
elseif strcmp(label,'14N')
    gamma=3.077E6;
elseif strcmp(label,'17O')
    gamma=-5.772E6;
elseif strcmp(label,'19F')
    gamma=40.078E6;
else
    error('unsupported nucleus label: %s',label);
end

end

function grumble(label)
if (~ischar(label))&&(~isstring(label))
    error('label must be a character string, or a string scalar.');
end
end

