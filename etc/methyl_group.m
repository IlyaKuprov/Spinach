% Coordinates for the four atoms of a methyl group. Syntax:
%
%         xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)
%
% Parameters:
%
%   c_xyz   - coordinates of C, row vector, Angstrom
%
%   cc_th   - polar theta angle of the C-C bond, radians
%
%   cc_ph   - polar phi angle of the C-C bond, radians
%
%   phase   - phase of the methyl group with respect to
%             its rotation around the C-C bond, radians
%
% Outputs:  
%
%     xyz   - a column cell array of Cartesian XYZ row
%             vectors; carbon is the first atom
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=methyl_group.m>

function xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)

% Check consistency
grumble(c_xyz,cc_th,cc_ph,phase);

% Generate a canonical methyl group
theta=acos(1/3);
xyz=[[0.000 0.000 0.000];
     [0.000 0.000 1.050]*euler2dcm(0,theta,0*pi/3+phase);
     [0.000 0.000 1.050]*euler2dcm(0,theta,2*pi/3+phase);
     [0.000 0.000 1.050]*euler2dcm(0,theta,4*pi/3+phase)];

% Rotate the CC bond
xyz=xyz*euler2dcm(0,cc_th,cc_ph);

% Translate the carbon
xyz=xyz+c_xyz;

% Return as a cell array
xyz=mat2cell(xyz,[1 1 1 1],3);

end

% Consistency enforcement
function grumble(c_xyz,cc_th,cc_ph,phase)
if (~isnumeric(c_xyz))||(~isreal(c_xyz))||(~isvector(c_xyz))||...
   (~isrow(c_xyz))||(numel(c_xyz)~=3)
    error('c_xyz mmust be a three-element real row vector.');
end
if (~isnumeric(cc_th))||(~isreal(cc_th))||(~isscalar(cc_th))
    error('c_th mmust be a real scalar.');
end
if (~isnumeric(cc_ph))||(~isreal(cc_ph))||(~isscalar(cc_ph))
    error('c_ph mmust be a real scalar.');
end
if (~isnumeric(phase))||(~isreal(phase))||(~isscalar(phase))
    error('phase mmust be a real scalar.');
end
end

% We've begun to long for the pitter-patter of little 
% feet - so we bought a dog. Well, it's cheaper, and
% you get more feet.
% 
% Rita Rudner

