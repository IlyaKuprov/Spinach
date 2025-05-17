% Returns true if none of the elements of the polyadic are Inf 
% or NaN. Syntax:
%
%                       answ=allfinite(p)
%
% Parameters:
%
%      p    - a polyadic object
%
% Outputs
%
%      answ - logical true if all numeric data
%             in the polyadic object is finite
% 
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/allfinite.m>

function answ=allfinite(p)

% Check the core array
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        if ~allfinite(p.cores{n}{k})
            answ=false; return
        end
    end
end

% Check prefix and suffix arrays
for n=1:numel(p.prefix)
    if ~allfinite(p.prefix{n})
        answ=false; return
    end
end
for n=1:numel(p.suffix)
    if ~allfinite(p.suffix{n})
        answ=false; return
    end
end

% All finite
answ=true();

end

% Beauty is the first test: there is no permanent 
% place in the world for ugly mathematics.
%
% G.H. Hardy

