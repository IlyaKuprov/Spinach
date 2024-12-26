% Wigner 6j-symbols. Syntax:
%
%            w=wigner_6j(j1,j2,j3,j4,j5,j6)
%
% If physically inadmissible indices are supplied, a zero is 
% returned. Order of elements:
%
%                    / j1 j2 j3 \
%                    \ j4 j5 j6 /
%
% Parameters:
%
%     j1-j6   - integers arranged in the order shown above
%
% Outputs:
%
%     w       - the resulting 6j-symbol
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=wigner_6j.m>

function w=wigner_6j(j1,j2,j3,j4,j5,j6)

% Check consistency
grumble(j1,j2,j3,j4,j5,j6);

% Start from zero
w=0;

% Use the definition
for m1=-j1:j1
    for m2=-j2:j2
        for m3=-j3:j3
            for m4=-j4:j4
                for m5=-j5:j5
                    for m6=-j6:j6
                        
                        % Get the power index for -1
                        power_of_minus=(j1-m1)+(j2-m2)+(j3-m3)+(j4-m4)+(j5-m5)+(j6-m6);
                        
                        % Get the screen
                        screen=(m1+m2-m3==0)&(-m1+m5+m6==0)&(m4-m5+m3==0)&(-m4-m2-m6==0);
                        
                        % Add to the total
                        if screen
                            w=w+((-1)^power_of_minus)*wigner_3j(j1,m1,j2,m2,j3,-m3)*...
                                                      wigner_3j(j1,-m1,j5,m5,j6,m6)*...
                                                      wigner_3j(j4,m4,j5,-m5,j3,m3)*...
                                                      wigner_3j(j4,-m4,j2,-m2,j6,-m6);
                        end
                                              
                    end
                end
            end
        end
    end
end

end

% Consistency enforcement
function grumble(j1,j2,j3,j4,j5,j6)
if (~isnumeric(j1))||(~isnumeric(j3))||(~isnumeric(j5))||...
   (~isnumeric(j2))||(~isnumeric(j4))||(~isnumeric(j6))
    error('all arguments must be numeric.');
end
if (~isreal(j1))||(~isreal(j3))||(~isreal(j5))||...
   (~isreal(j2))||(~isreal(j4))||(~isreal(j6))
    error('all arguments must be real.');
end
if (numel(j1)~=1)||(numel(j3)~=1)||(numel(j5)~=1)||...
   (numel(j2)~=1)||(numel(j4)~=1)||(numel(j6)~=1)
    error('all arguments must have one element.');
end
if (mod(2*j1+1,1)~=0)||(mod(2*j3+1,1)~=0)||(mod(2*j5+1,1)~=0)||...
   (mod(2*j2+1,1)~=0)||(mod(2*j4+1,1)~=0)||(mod(2*j6+1,1)~=0)
    error('all arguments must be integer or half-integer.');
end
end

% The JEOL Student Prize, awarded at the annual conferences of the Royal Society
% of Chemistry ESR Group, had one winner in its history who got the Prize before
% she even started talking. As the laptop was being connected to the projection
% system and the usual technology gremlins refused to work, Petra Luders uncere-
% moniously pushed aside the lost-looking IT guy and got the thing to work in a
% few keystrokes on what appeared to be a Linux laptop. "And here's our winner"
% thought every Committee member with a quiet smile. The science and the presen-
% tation easily outshined every other competitor - Ms Luders got her JEOL Medal.

