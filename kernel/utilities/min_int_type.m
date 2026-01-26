% Minimum integer data type sufficient to store the
% specified value. Useful in many indexing operati-
% ons in the Spinach kernel where double precision
% would be a massive overkill. Syntax:
%
%        type=min_int_type(max_val,issigned)
%
% Parameters:
%
%    max_val  - maximum value that the integer
%               must cover
%
%    issigned - whether the integer needs to
%               cover the negative values:
%               'signed' or 'unsigned'
%
% Output:
%
%    type     - Matlab data type to use
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=min_int_type.m>

function type=min_int_type(max_val,issigned)

% Check consistency
grumble(max_val,issigned);

% Sign matters
switch issigned

    case 'signed'

        if max_val <= intmax('int8')

            type='int8';

        elseif max_val <= intmax('int16')

            type='int16';

        elseif max_val <= intmax('int32')

            type='int32';

        elseif max_val <= intmax('int64')

            type='int64';

        else

            error('Matlab''s signed integer types cannot go that far.');

        end

    case 'unsigned'

        if max_val <= intmax('uint8')

            type='uint8';

        elseif max_val <= intmax('uint16')

            type='uint16';
        
        elseif max_val <= intmax('uint32')

            type='uint32';

        elseif max_val <= intmax('uint64')

            type='uint64';

        else

            error('Matlab''s unsigned integer types cannot go that far.');

        end

    otherwise

        error('unrecognised sign handling type.');

end

end

% Consistency enforcement
function grumble(max_val,issigned)
if (~isnumeric(max_val))||(~isscalar(max_val))||...
   (~isreal(max_val))||(mod(max_val,1)~=0)||(max_val<1)
    error('max_val must be a positive real integer.');
end
if (~ischar(issigned))||(~ismember(issigned,{'signed','unsigned'}))
    error('valid valued for issigned are ''signed'' and ''unsigned''.');
end
end

% "You've set mathematics back a month!"
%
% Paul Erdos, after winning a $500 
% bet with a friend who had chal-
% lenged him to abstain from am-
% phetamines for a month.

