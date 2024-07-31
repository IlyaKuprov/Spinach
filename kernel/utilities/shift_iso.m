% Replaces the isotropic parts of interaction tensors with user-
% supplied values. This is useful for correcting DFT calculations,
% where the anisotropy of the various spin interactions is usual-
% ly satisfactory, but the isotropic part is not. Syntax:
%
%        tensors=shift_iso(tensors,spin_numbers,new_iso)
%
% Parameters:
%
%      tensors      - a cell array of interaction tensors 
%                     as 3x3 matrices
%
%      spin_numbers - a vector containing the numbers 
%                     of spins in the tensors array that
%                     should have the isotropic parts
%                     replaced
%
%      new_iso      - a vector containing the new isotro-
%                     pic parts in the same order as the
%                     spin numbers listed in spin_numbers
%
% Outputs:
%
%      tensors      - a cell array of interaction tensors 
%                     as 3x3 matrices 
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=shift_iso.m>

function tensors=shift_iso(tensors,spin_numbers,new_iso)

% Check consistency
grumble(tensors,spin_numbers,new_iso)

% Loop over the tensors
for n=1:numel(spin_numbers)
    
    % Isolate the anisotropy
    [~,rank1,rank2]=mat2sphten(tensors{spin_numbers(n)});
    
    % Rebuild with the new isotropic part
    tensors{spin_numbers(n)}=sphten2mat([],rank1,rank2)+new_iso(n)*eye(3);
    
end

end

% Consistency enforcement
function grumble(tensors,spin_numbers,new_iso)
if (~iscell(tensors))||any(any(~cellfun(@isreal,tensors)))||...
   any(any(~cellfun(@(x)all(size(x)==[3 3]|isempty(x)),tensors)))
    error('tensors parameter must be a cell array of real 3x3 matrices.');
end
if (~isnumeric(spin_numbers))||any(mod(spin_numbers,1)~=0)||any(spin_numbers<1)
    error('spin_numbers must be a vector of positive integers.');
end
if any(spin_numbers>numel(tensors))
    error('an index in spin_numbers is greater than the number of tensors supplied.');
end
if (~isnumeric(new_iso))||(~isreal(new_iso))
    error('new_iso must be a vector of real numbers.');
end
if numel(spin_numbers)~=numel(new_iso)
    error('the number of elements in spin_number and new_iso parameters must be the same.');
end
end

% There once was an X from place B,
% Who satisfied predicate P,
% Then X did thing A,
% In a specified way,
% Resulting in circumstance C.

