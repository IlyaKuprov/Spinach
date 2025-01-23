% Estimates the rotational correlation time from the 
% longitudinal CSA relaxation rate. Syntax:
%
%        tauc=r1csa2tauc(R1,del_sq,B0,isotope)
%
% Parameters:
%
%    R1      - longitudinal relaxation rate, Hz
%
%    del_sq  - second rank invariant of the CSA,
%              see blinv.m function
%             
%    B0      - magnetic field, Tesla
%
%    isotope - isotope specification string, e.g. '1H'
%
% Outputs:
%
%    tauc    - rotational correlation time, seconds
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=r1csa2tauc.m>

function tauc=r1csa2tauc(R1,del_sq,B0,isotope)

% Check consistency
grumble(R1,del_sq,B0,isotope);

% Get the Zeeman frequency
omega=-B0*spin(isotope);

% Solve the quadratic equation
tauc(1)=(del_sq*omega^2-sqrt(del_sq^2*omega^4-225*omega^2*R1^2))/(15*omega^2*R1);
tauc(2)=(del_sq*omega^2+sqrt(del_sq^2*omega^4-225*omega^2*R1^2))/(15*omega^2*R1);

% Check for imaginary components
if (~isreal(tauc(1)))&&(~isreal(tauc(2)))
    error('no real solutions - physically impossible case.');
end

end

% Consistency enforcement
function grumble(R1,del_sq,B0,isotope)
if (~isnumeric(R1))||(~isreal(R1))||...
   (~isscalar(R1))||(R1<=0)
    error('R1 must be a positive real number.');
end
if (~isnumeric(del_sq))||(~isreal(del_sq))||...
   (~isscalar(del_sq))||(del_sq<=0)
    error('del_sq must be a positive real number.');
end
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))
    error('B0 must be a real number.');
end
if ~ischar(isotope)
    error('isotope must be a character string.');
end
end

% Yet again, the attack is described as 'cowardly'. This 
% is simply untrue: it must require immense, though repel-
% lent, courage to blow yourself up.
%
% Charles Moore

