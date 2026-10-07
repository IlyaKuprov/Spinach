% Leakage-based time estimate for zero track elimination. Syntax:
%
%          duration=zte_warr(spin_system,L,rho,projector)
%
% Parameters:
%
%    spin_system - Spinach structure with the positive absolute
%                  state 2-norm tolerance tols.zte_warr
%
%    L           - fixed square Liouvillian, in radians per second
%
%    rho         - initial state column vector
%
%    projector   - coordinate embedding returned by zte, or the
%                  scalar 1 when no reduction is performed
%
% Output:
%
%    duration    - estimated time in seconds before the absolute
%                  state 2-norm error reaches tols.zte_warr
%
% Notes: for retained coordinates S and discarded coordinates D,
%        the estimate is max(0,tolerance-d)/(b*r), with d=norm(rho(D)),
%        r=norm(rho(S)), and b=norm(L(D,S),'fro'). The Frobenius norm
%        cheaply bounds the leakage block's spectral norm. This
%        linear leakage estimate neglects subsequent amplification;
%        it is not guaranteed under amplifying dynamics.
%        Initial discard at or above tolerance gives zero seconds;
%        zero estimated leakage otherwise gives Inf. No reduction
%        or a zero input state gives Inf. Only the supplied vector
%        and fixed generator are covered, not a slowpass spectrum,
%        numerical propagation error, or roundoff.
%
% ilya.kuprov@weizmann.ac.il

function duration=zte_warr(spin_system,L,rho,projector)

% Validate the input
grumble(spin_system,L,rho,projector);

% Identify the retained coordinates without expanding a scalar identity
if isequal(projector,1)
    retained=true(size(rho));
else
    retained=any(projector,2);
end

% Estimate the initial error and the retained-to-discarded leakage rate
discarded=norm(rho(~retained));
rate=norm(L(~retained,retained),'fro')*norm(rho(retained));
duration=Inf;
if discarded>=spin_system.tols.zte_warr
    duration=0;
elseif rate>0
    duration=(spin_system.tols.zte_warr-discarded)/rate;
end

% Report the estimate without claiming a error guarantee
report(spin_system,['ZTE warranty estimate: absolute 2-norm tolerance ' ...
                    num2str(spin_system.tols.zte_warr,17) ', estimated time ' ...
                    num2str(duration,17) ' seconds (supplied vector, fixed generator; not a guarantee).']);

end

% Input validation function
function grumble(spin_system,L,rho,projector)
if (~isnumeric(L))||(~ismatrix(L))||(size(L,1)~=size(L,2))||...
   any(~isfinite(nonzeros(L)))
    error('L must be a finite square numeric matrix.');
end
if (~isnumeric(rho))||(size(rho,2)~=1)||(size(rho,1)~=size(L,1))||...
   any(~isfinite(nonzeros(rho)))
    error('rho must be a finite column vector matching L.');
end
if (~isnumeric(projector))||(~ismatrix(projector))||...
   (~isequal(projector,1)&&(size(projector,1)~=size(L,1)))
    error('projector must be the coordinate embedding returned by zte.');
end
if (~isnumeric(spin_system.tols.zte_warr))||(~isreal(spin_system.tols.zte_warr))||...
   (~isscalar(spin_system.tols.zte_warr))||(~isfinite(spin_system.tols.zte_warr))||...
   (spin_system.tols.zte_warr<=0)
    error('zte_warr must be a finite positive real scalar.');
end
end


