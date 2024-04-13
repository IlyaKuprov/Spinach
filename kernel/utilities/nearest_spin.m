% Returns the index of the nearest spin to the one speci-
% fied. Only spins for which Cartesian coordinates are
% available are considered. Syntax:
%
%              k=nearest_spin(spin_system,n)
%
% Parameters:
%
%    n   - index of the spin in question
%
% Outputs:
%
%    k   - index of the nearest spin
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=nearest_spin.m>

function [k,d]=nearest_spin(spin_system,n)

% Check consistency
grumble(spin_system,n);

% Starting point
k=[]; d=inf;

% Find the nearest spin
for s=1:numel(spin_system.inter.coordinates)
    if (s~=n)&&(~isempty(spin_system.inter.coordinates{s}))
        current_dist=norm(spin_system.inter.coordinates{s}-... % #NORMOK
                          spin_system.inter.coordinates{n},2);
        if current_dist<d, k=s; d=current_dist; end
    end
end

% Catch pathological cases
if isempty(k), error('no other spin has coordinates'); end

end

% Consistency enforcement
function grumble(spin_system,n)
if (~isnumeric(n))||(~isscalar(n))||...
   (~isreal(n))||(n<1)||(mod(n,1)~=0)
    error('n must be a positive real integer');
end
if n>numel(spin_system.comp.isotopes)
    error('the specified spin does not exist');
end
if isempty(spin_system.inter.coordinates{n})
    error('the specified spin does not have coordinates');
end
end

% К этому времени мои антисоветские стихи приумножились;
% таланта в них не прибавилось, но как листовки они смо-
% трелись. Набирая свое тайное общество, я всем встреч-
% ным и поперечным их давала читать. В ИНЯЗе работали и
% учились редкие люди: опять никто не донёс!
%
% Валерия Новодворская

