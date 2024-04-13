% Linear least squares fitter for residual dipolar couplings. Isotope 
% pairs are arbitrary; multiple isotope pairs may be supplied at the
% same time. Syntax:
%
%                     S=rdc_fit(isotopes,xyz,rdc)
%
% Parameters:
%
%   isotopes - N x 2 cell array of strings with Spinach 
%              isotope specifications, e.g. '13C'
%
%   xyz      - N x 2 cell array of 3-element vectors
%              with Cartesian coordinates in Angstrom
%
%   rdc      - N x 1 vector with residual dipolar coup-
%              lings in Hz
%
% Outputs:
%
%   S        - Saupe order matrix, a symmetric real
%              dimensionless 3x3 matrix
%                
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rdc_fit.m>

function S=rdc_fit(isotopes,xyz,rdc)

% Check consistency
grumble(isotopes,xyz,rdc);

% Build dipolar tensor array
D=cell(1,size(isotopes,1));
for n=1:size(isotopes,1)

    % Get the dipole-dipole coupling tensor
    [~,~,~,~,D{n}]=xyz2dd(xyz{n,1},xyz{n,2},isotopes{n,1},isotopes{n,2});

    % Extract second rank component
    [~,~,D{n}]=mat2sphten(D{n});

end

% Run pseudoinverse least squares
S=(3/2)*(cell2mat(D)'\(2*pi*rdc(:)));

% Convert back into a matrix
S=real(sphten2mat([],[],S));

end

% Consistency enforcement
function grumble(isotopes,xyz,rdc)
if (~iscell(isotopes))||(size(isotopes,2)~=2)
    error('isotopes must be an N x 2 cell array of character strings.');
end
for n=1:size(isotopes,1)
    if (~ischar(isotopes{n,1}))||(~ischar(isotopes{n,2}))
        error('isotopes must be an N x 2 cell array of character strings.');
    end
    if strcmp(isotopes{n,1},isotopes{n,2})
        error('spins in each pair must belong to different isotopes.');
    end
end
if (~iscell(xyz))||(size(xyz,2)~=2)
    error('xyz must be an N x 2 cell array of three-element vectors.');
end
for n=1:size(xyz,1)
    if (~isnumeric(xyz{n,1}))||(~isreal(xyz{n,1}))||(numel(xyz{n,1})~=3)||...
       (~isnumeric(xyz{n,2}))||(~isreal(xyz{n,2}))||(numel(xyz{n,2})~=3)
        error('elements of xyz must be real three-element vectors');
    end
end
if (~isnumeric(rdc))||(~isreal(rdc))||(~iscolumn(rdc))
    error('rdc must be a column vector with real elements.');
end
if (size(isotopes,1)~=size(xyz,1))||(size(xyz,1)~=size(rdc,1))
    error('all inputs must have the same number of rows.');
end
end

% Эт конечно, царь мне вставит, у него эт не отнять. 
% Но однако, слышь ты зелень, как бы мне тебе сказать... 
% Лучше получу я в ухо, лучше в челюсть получу,
% Но зато я буду холост, буду жить я как хочу.
% 
% Сектор Газа, "Ария Ивана и Лягушки"

