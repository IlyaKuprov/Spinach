%SPIN_LABEL Converts Kehl isotope labels into Spinach isotope labels.
%
%   LABEL=SPIN_LABEL(LABEL_IN) returns a character vector suitable for
%   Spinach spin(), create(), and basis() calls. The only legacy alias
%   currently needed is 2D -> 2H.


function label=kehl_spin_label(label_in)

    % Check consistency
    grumble(label_in);
if isstring(label_in)
    label=char(label_in);
else
    label=label_in;
end

if strcmp(label,'2D')
    label='2H';
end

end

function grumble(label_in)
if (~ischar(label_in))&&(~isstring(label_in))
    error('label_in must be a character string, or a string scalar.');
end
end

