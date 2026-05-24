% Real roots of a cubic polynomial in the unit interval. Syntax:
%
%           root_list=cubic_roots(poly_coeffs,root_tol)
%
% Parameters:
%
%     poly_coeffs - four real coefficients [a b c d] of
%                   a*x^3+b*x^2+c*x+d
%
%        root_tol - positive real root filtering tolerance
%
% Outputs:
%
%     root_list   - sorted row vector of real roots in [0,1]
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cubic_roots.m>

function root_list=cubic_roots(poly_coeffs,root_tol)

% Check consistency
grumble(poly_coeffs,root_tol);

% Normalise polynomial coefficients
poly_coeffs=poly_coeffs(:).';
poly_scale=max(abs(poly_coeffs));
if poly_scale==0
    root_list=[];
    return;
end
poly_coeffs=poly_coeffs/poly_scale;

% Drop leading numerical zeros
lead_idx=find(abs(poly_coeffs)>root_tol,1,'first');
if isempty(lead_idx)
    root_list=[];
    return;
end
poly_coeffs=poly_coeffs(lead_idx:end);

% Find real roots inside the unit interval
root_list=roots(poly_coeffs);
root_list=root_list(abs(imag(root_list))<root_tol);
root_list=real(root_list);
root_list=root_list((root_list>=-root_tol)&...
                    (root_list<=1+root_tol));
root_list=min(max(root_list,0),1);
root_list=sort(root_list(:).');

% Merge numerically coincident roots
if ~isempty(root_list)
    root_list=root_list([true abs(diff(root_list))>root_tol]);
end

end

% Consistency enforcement
function grumble(poly_coeffs,root_tol)
if (~isnumeric(poly_coeffs))||(~isreal(poly_coeffs))||...
   (numel(poly_coeffs)~=4)||any(~isfinite(poly_coeffs(:)))
    error('poly_coeffs must be a finite real vector with four elements.');
end
if (~isnumeric(root_tol))||(~isreal(root_tol))||...
   (~isscalar(root_tol))||(~isfinite(root_tol))||(root_tol<=0)
    error('root_tol must be a finite positive real scalar.');
end
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

