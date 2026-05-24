% Real roots of a cubic polynomial in the unit interval. Syntax:
%
%                    root_list=cubic_roots(poly_coeffs,root_tol)
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


