% Returns true if the polyadic representation is real. Syntax:
%
%                         answer=isreal(p)
%
% Parameters:
%
%     p      - a polyadic object
%
% Outputs:
%
%     answer - true if all numeric data in the object is real
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/isreal.m>

function answer=isreal(p)

% Check the core array
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        if ~isreal(p.cores{n}{k})
            answer=false; return
        end
    end
end

% Check prefix and suffix arrays
for n=1:numel(p.prefix)
    if ~isreal(p.prefix{n})
        answer=false; return
    end
end
for n=1:numel(p.suffix)
    if ~isreal(p.suffix{n})
        answer=false; return
    end
end

% All data is real
answer=true();

end

% A little inaccuracy sometimes saves a ton of explanation.
%
% H.H. Munro

