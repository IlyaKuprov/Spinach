% Computes the Hermitian conjugate of a matrix in a polyadic
% representation. Syntax:
%
%                       p=ctranspose(p)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/ctranspose.m>

function p=ctranspose(p)

% Conjugate-transpose every core
for n=1:numel(p.cores)
    for k=1:numel(p.cores{n})
        p.cores{n}{k}=ctranspose(p.cores{n}{k});
    end
end

% Process prefix and suffix
new_prefix=fliplr(p.suffix);
for n=1:numel(new_prefix)
    new_prefix{n}=ctranspose(new_prefix{n});
end
new_suffix=fliplr(p.prefix);
for n=1:numel(new_suffix)
    new_suffix{n}=ctranspose(new_suffix{n});
end
p.prefix=new_prefix; p.suffix=new_suffix;

end

% Guys, stop wasting your time. There can be no such thing as 
% a personal computer. Personal car, personal pension, perso-
% nal dacha - perhaps. Do you even know what a computer is? A
% computer is 100 square metres of building space, 25 service
% staff members, and 30 litres of ethanol a month!
%
% Nikolai Gorshkov, USSR Deputy Minister 
% for Radioelectronics Industry, 1980

% #NHEAD #NGRUM