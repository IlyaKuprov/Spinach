% Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri-
% cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2)
% spherical harmonics. In the more general sense, the coefficient refers
% to the expansion coefficient of |L,M> angular momentum or spin state
% in the product basis of |L1,M1>|L2,M2> states. Syntax:
%
%                       cg=cg_fast(L,M,L1,M1,L2,M2)
%
% Parameters:
%
%     L,M,L1,M1,L2,M2  - integer or half-integer indices of 
%                        the angular momentum or spin states
%
% Outputs:
%
%     cg               - floating-point (double precision) 
%                        Clebsch-Gordan coefficient
%                                
% Note: only some combinations of L,M,L1,M1,L2,M2 are allowed by the pro-
%       perties of spherical harmonics and spin states. If inadmissible
%       indices are supplied, zero is returned.
%
% Note: CG coefficient calculation in double-precision arithmetic is not
%       a trivial matter for high ranks. This function produces fast ans-
%       wers with an accuracy of about 1e-3 up to about L=20. A slower 
%       machine precision implementation for higher ranks is available
%       in clebsch_gordan.m function.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cg_fast.m>

function cg=cg_fast(L,M,L1,M1,L2,M2)

% Check consistency
grumble(L,M,L1,M1,L2,M2);

% Match the notation to Varshalovich, Section 8.2.1
cg=0; c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2;

% Run zero tests (Stage I)
prefactor_is_nonzero=(a+alp>=0)&&(a-alp>=0)&&...
                     (b+bet>=0)&&(b-bet>=0)&&...
                     (c+gam>=0)&&(c-gam>=0)&&...
                     (gam==alp+bet);
    
% Proceed if appropriate
if prefactor_is_nonzero
    
    % Run zero tests (Stage II)
    delta_is_nonzero=(a+b-c>=0)&&(a-b+c>=0)&&...
                     (-a+b+c>=0)&&(a+b+c+1>=0);
                     
    % Proceed if appropriate
    if delta_is_nonzero
            
        % Run zero tests (Stage III)
        lower_sum_limit=max([alp-a, b+gam-a, 0]);
        upper_sum_limit=min([c+b+alp, c+b-a, c+gam]);
        sum_is_nonzero=(upper_sum_limit>=lower_sum_limit);
    
        % Proceed if appropriate
        if sum_is_nonzero
            
            % Compute Equation 8.2.1(5) using log-factorials
            for z=lower_sum_limit:upper_sum_limit
                cg=cg+((-1)^(b+bet+z))*sqrt(2*c+1)*...
                   exp(logfactorial(c+b+alp-z)+logfactorial(a-alp+z)-logfactorial(z)-...
                       logfactorial(c-a+b-z)-logfactorial(c+gam-z)-logfactorial(a-b-gam+z)+...
                      (logfactorial(a+b-c)+logfactorial(a-b+c)+logfactorial(-a+b+c)+...
                       logfactorial(c+gam)+logfactorial(c-gam)-logfactorial(a+b+c+1)-...
                       logfactorial(a+alp)-logfactorial(a-alp)-logfactorial(b+bet)-...
                       logfactorial(b-bet))/2);
            end
            
        end
        
    end
    
end

end

% Consistency enforcement
function grumble(L,M,L1,M1,L2,M2)
if (~isnumeric(L))||(~isnumeric(L1))||(~isnumeric(L2))||...
   (~isnumeric(M))||(~isnumeric(M1))||(~isnumeric(M2))
    error('all arguments must be numeric.');
end
if (~isreal(L))||(~isreal(L1))||(~isreal(L2))||...
   (~isreal(M))||(~isreal(M1))||(~isreal(M2))
    error('all arguments must be real.');
end
if (numel(L)~=1)||(numel(L1)~=1)||(numel(L2)~=1)||...
   (numel(M)~=1)||(numel(M1)~=1)||(numel(M2)~=1)
    error('all arguments must have one element.');
end
if (mod(2*L+1,1)~=0)||(mod(2*L1+1,1)~=0)||(mod(2*L2+1,1)~=0)||...
   (mod(2*M+1,1)~=0)||(mod(2*M1+1,1)~=0)||(mod(2*M2+1,1)~=0)
    error('all arguments must be integer or half-integer.');
end
if (L>20)||(L1>20)||(L2>20)
    error('the answer will not be accurate, use clebsch_gordan() instead.');
end
end

% If you follow the trodden path, you may find
% that all the grass has been eaten.
%
% Andre Geim

