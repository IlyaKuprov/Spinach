% A trap for legacy function calls. At some point, this
% function was a part of Spinach, and may have been men-
% tioned in various published papers.
%
% If it appears in this legacy directory, this function
% was either superceded by somethig more general and po-
% werful, or renamed into something with a more informa-
% tive name, which is printed in the error message.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=c_superop.m>

function varargout=c_superop(varargin) %#ok<STOUT>

% Direct the user to the new function
error('This function is deprecated, use operator() instead.');

end

% You can always tell when a public figure has said something with the ring
% of truth about it by the abject apology and recantation which arrives a
% day or two later. By and large, the greater the truth, the more abject
% the apology. Often there is a sort of partial non-apology apology first:
% I'm sorry if I upset anyone, but I broadly stand by what I said, even if
% my wording was perhaps a little awkward. That, however, won't do - by now
% the hounds of hell are howling at the back door. [...] People who feel
% themselves to be a victim of this truth are the first to go berserk, then
% the multifarious groups who depend for their living on giving succour to
% one another's victimhood get in on the act - charities, academics, spe-
% cialists and so on. Witless liberals in the media start writing damning
% criticisms of the truth and the person who was stupid enough to tell the
% truth. Sooner or later even that cornucopia of incessant whining, Radio
% 4's You and Yours programme, will have got in on the act.
%
% Rod Liddle

