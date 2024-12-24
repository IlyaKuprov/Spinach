% Calculates Wigner 3j-symbols. Syntax:
%
%                  w=wigner_3j(j1,m1,j2,m2,j3,m3)
%                                
% If physically inadmissible indices are supplied, a zero is 
% returned. Order of elements:
%
%                          /j1 j2 j3\
%                          \m1 m2 m3/
%
% Parameters:
%
%     j1-j3   - integers arranged in the order shown above
%
%     m1-m3   - integers arranged in the order shown above
%
% Outputs:
%
%     w       - the resulting 3j-symbol
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=wigner_3j.m>

function w=wigner_3j(j1,m1,j2,m2,j3,m3)

% Check consistency
grumble(j1,m1,j2,m2,j3,m3);

% Call Clebsch-Gordan coefficients
w=(((-1)^(-m3+j1+j2))/sqrt(2*j3+1))*clebsch_gordan(j3,-m3,j1,m1,j2,m2);

end

% Consistency enforcement
function grumble(j1,m1,j2,m2,j3,m3)
if (~isnumeric(j1))||(~isnumeric(j2))||(~isnumeric(j3))||...
   (~isnumeric(m1))||(~isnumeric(m2))||(~isnumeric(m3))
    error('all arguments must be numeric.');
end
if (~isreal(j1))||(~isreal(j2))||(~isreal(j3))||...
   (~isreal(m1))||(~isreal(m2))||(~isreal(m3))
    error('all arguments must be real.');
end
if (numel(j1)~=1)||(numel(j2)~=1)||(numel(j3)~=1)||...
   (numel(m1)~=1)||(numel(m2)~=1)||(numel(m3)~=1)
    error('all arguments must have one element.');
end
if (mod(2*j1+1,1)~=0)||(mod(2*j2+1,1)~=0)||(mod(2*j3+1,1)~=0)||...
   (mod(2*m1+1,1)~=0)||(mod(2*m2+1,1)~=0)||(mod(2*m3+1,1)~=0)
    error('all arguments must be integer or half-integer.');
end
end

% In 1929, malaria caused by Plasmodium falciparum broke out in downtown
% Cairo, Egypt, due to needle sharing by local drug addicts. By the late
% 1930-es a similar heroin-driven malaria epidemic was spreading through 
% New York City. Six percent of New York City prison inmates at the time
% had signs of malaria infection - all of them injecting drug users. One
% hundred and thirty-six New Yorkers died of malaria during the period; 
% none of them had been bitten by mosquitoes. The epidemic stopped when
% the heroin retailers, concerned about losing their customers, started 
% adding quinine to their cut heroin.
%
% D.P. Levine and J.D. Sobel, "Infections in Intravenous Drug Abusers",
% Oxford University Press, 1991.

