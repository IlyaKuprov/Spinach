% Spin system of strychnine. Isotropic chemical shifts and J-couplings
% are taken from "200 and more NMR experiments: a practical course" by
% by Berger and Braun. Coordinates taken as those of the major confor-
% mer of strychnine proposed in http://dx.doi.org/10.1039/C0CC04114A
% Syntax:
%
%                   [sys,inter]=strychnine(spins)
%
% The 'spins' parameter is a cell array containing the isotopes to im-
% port, e.g. {'1H','13C'}.
%
% Note: 13C-13C J-couplings are not provided - this file is only
%       suitable for natural abundance 13C simulations.
%
% Note: CSA tensors are not provided - relaxation theory treatments
%       on top of this file would not account for CSA relaxation.
%
% Note: only shifts and coordinates are provided for 15N nuclei.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=strychnine.m>

function [sys,inter]=strychnine(spins)

% Check consistency
grumble(spins);

% Shorthands for human-readable coupling designations below
H1=1; H2=2; H3=3; H4=4; H8=5; H11a=6; H11b=7; H12=8; H13=9; H14=10;
H15a=11; H15b=12; H16=13; H17a=14; H17b=15; H18a=16; H18b=17; H20a=18;
H20b=19; H22=20; H23a=21; H23b=22;

% Hydrogen atoms
H_isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
            '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};
numH=numel(H_isotopes);

% Proton chemical shifts
H_zeeman=cell(1,numH);
H_zeeman{H3}=  7.255;  H_zeeman{H22}= 5.915;
H_zeeman{H2}=  7.098;  H_zeeman{H1}=  7.167;
H_zeeman{H4}=  8.092;  H_zeeman{H12}= 4.288;
H_zeeman{H23b}=4.066;  H_zeeman{H23a}=4.148;
H_zeeman{H16}= 3.963;  H_zeeman{H8}=  3.860;
H_zeeman{H20b}=2.745;  H_zeeman{H20a}=3.716;
H_zeeman{H18a}=3.219;  H_zeeman{H18b}=2.878;
H_zeeman{H13}= 1.276;  H_zeeman{H17a}=1.890;
H_zeeman{H17b}=1.890;  H_zeeman{H11a}=3.132;
H_zeeman{H11b}=2.670;  H_zeeman{H14}= 3.150;
H_zeeman{H15b}=1.462;  H_zeeman{H15a}=2.360;

% Proton J-couplings
H_coupling=cell(numH);
H_coupling{H3,H4}=      7.90;  H_coupling{H22,H23a}=   7.00;
H_coupling{H22,H23b}=   6.10;  H_coupling{H2,H3}=      7.44;
H_coupling{H2,H4}=      0.98;  H_coupling{H1,H2}=      7.49;
H_coupling{H1,H3}=      1.08;  H_coupling{H1,H4}=      0.23;
H_coupling{H12,H13}=    3.30;  H_coupling{H23a,H23b}=-13.80;
H_coupling{H8,H13}=    10.41;  H_coupling{H20a,H20b}=-14.80;
H_coupling{H20a,H22}=   1.79;  H_coupling{H18a,H18b}=-13.90;
H_coupling{H13,H14}=    3.29;  H_coupling{H17a,H17b}=-13.90;
H_coupling{H17a,H18a}=  5.50;  H_coupling{H17a,H18b}=  7.20;
H_coupling{H17b,H18a}=  3.20;  H_coupling{H17b,H18b}= 10.70;
H_coupling{H11a,H11b}=-17.34;  H_coupling{H11a,H12}=   3.34;
H_coupling{H11b,H12}=   8.47;  H_coupling{H14,H15a}=   4.11;
H_coupling{H14,H15b}=   1.96;  H_coupling{H14,H22}=    0.47;
H_coupling{H14,H20a}=   1.61;  H_coupling{H15a,H15b}=-14.35;
H_coupling{H15a,H16}=   4.33;  H_coupling{H15b,H16}=   2.42;

% Proton coordinates
H_coordinates=cell(1,numH);
H_coordinates{H1}=  [ 2.7336, -2.8961, -0.4256];
H_coordinates{H2}=  [ 5.0747, -2.1373, -0.7803];
H_coordinates{H3}=  [ 5.6630,  0.2547, -0.5015];
H_coordinates{H4}=  [ 3.9144,  1.9280,  0.1137];
H_coordinates{H8}=  [-0.6004,  0.4789,  1.3185];
H_coordinates{H11a}=[-0.3262,  4.0926, -0.0781];
H_coordinates{H11b}=[-0.6863,  3.1492,  1.3616];
H_coordinates{H12}= [-1.7489,  2.7068, -1.3360];
H_coordinates{H13}= [-0.0861,  1.0611, -1.5871];
H_coordinates{H14}= [-2.2234,  0.1869, -2.2835];
H_coordinates{H15a}=[ 0.0099, -1.0345, -2.4296];
H_coordinates{H15b}=[-1.3580, -2.1373, -2.4297];
H_coordinates{H16}= [ 0.3886, -2.9088, -0.8537];
H_coordinates{H17a}=[ 1.0318, -2.7083,  1.7742];
H_coordinates{H17b}=[ 0.6059, -1.1921,  2.5799];
H_coordinates{H18a}=[-1.2874, -3.1432,  2.3544];
H_coordinates{H18b}=[-1.7510, -1.4610,  2.0500];
H_coordinates{H20a}=[-2.9738, -2.8430, -1.0448];
H_coordinates{H20b}=[-3.4330, -2.5457,  0.6228];
H_coordinates{H22}= [-4.4597, -0.3928,  0.8115];
H_coordinates{H23a}=[-3.8532,  1.8278, -1.0768];
H_coordinates{H23b}=[-4.5333,  1.9818,  0.5498];

% Carbon and nitrogen atoms
CN_isotopes={'13C','13C','13C','13C','13C','13C','13C','13C','15N','13C',...
             '13C','13C','13C','13C','13C','13C','13C','13C','15N','13C',...
             '13C','13C','13C'};
numCN=numel(CN_isotopes);

% Carbon and nitrogen chemical shifts
CN_zeeman=cell(1,numCN);
CN_zeeman{10}= 169.28; CN_zeeman{9}=  155.00;
CN_zeeman{5}=  142.23; CN_zeeman{21}= 140.45;
CN_zeeman{6}=  132.72; CN_zeeman{3}=  128.56;
CN_zeeman{22}= 127.34; CN_zeeman{2}=  124.20;
CN_zeeman{1}=  122.26; CN_zeeman{4}=  116.23;
CN_zeeman{12}=  76.85; CN_zeeman{23}=  64.60;
CN_zeeman{16}=  60.28; CN_zeeman{8}=   60.10;
CN_zeeman{20}=  52.68; CN_zeeman{7}=   51.96;
CN_zeeman{18}=  50.35; CN_zeeman{13}=  48.22;
CN_zeeman{17}=  42.85; CN_zeeman{11}=  42.48;
CN_zeeman{19}=  35.00; CN_zeeman{14}=  31.60;
CN_zeeman{15}=  26.84;

% Carbon and nitrogen coordinates
CN_coordinates=cell(1,numCN+2);
CN_coordinates{4}= [ 3.6711,  0.8813, -0.0049];
CN_coordinates{3}= [ 4.6363, -0.0692, -0.3529];
CN_coordinates{2}= [ 4.3069, -1.4174, -0.5113];
CN_coordinates{1}= [ 2.9885, -1.8448, -0.3155];
CN_coordinates{6}= [ 2.0128, -0.9131,  0.0231];
CN_coordinates{5}= [ 2.3591,  0.4385,  0.1719];
CN_coordinates{8}= [ 0.0012,  0.3360,  0.4149];
CN_coordinates{7}= [ 0.5532, -1.1283,  0.3650];
CN_coordinates{9}= [ 1.2080,  1.1956,  0.4953];
CN_coordinates{17}=[ 0.3655, -1.8391,  1.7302];
CN_coordinates{18}=[-1.0998, -2.2803,  1.7028];
CN_coordinates{19}=[-1.3733, -2.6144,  0.2962];
CN_coordinates{16}=[-0.2694, -2.0788, -0.5579];
CN_coordinates{15}=[-0.8203, -1.4052, -1.8175];
CN_coordinates{14}=[-1.7397, -0.2431, -1.3932];
CN_coordinates{13}=[-0.8137,  0.8361, -0.7930];
CN_coordinates{10}=[ 1.0821,  2.5627,  0.3141];
CN_coordinates{11}=[-0.3587,  3.0755,  0.3170];
CN_coordinates{12}=[-1.4198,  2.2035, -0.4123];
CN_coordinates{23}=[-3.7138,  1.5364, -0.0219];
CN_coordinates{22}=[-3.7170,  0.0350,  0.1391];
CN_coordinates{21}=[-2.8213, -0.7674, -0.4478];
CN_coordinates{20}=[-2.7216, -2.2506, -0.1548];

% Combine the arrays
sys.isotopes=[H_isotopes, CN_isotopes];
inter.zeeman.scalar=[H_zeeman, CN_zeeman];
inter.coordinates=[H_coordinates, CN_coordinates];
inter.coupling.scalar=cell(numH+numCN);
inter.coupling.scalar(1:numH,1:numH)=H_coupling;

% Add 13C-1H J-couplings
inter.coupling.scalar{H3,3+numH}=   159.6;
inter.coupling.scalar{H22,22+numH}= 157.7;
inter.coupling.scalar{H2,2+numH}=   160.9;
inter.coupling.scalar{H1,1+numH}=   159.0;
inter.coupling.scalar{H4,4+numH}=   168.0;
inter.coupling.scalar{H12,12+numH}= 145.4;
inter.coupling.scalar{H23b,23+numH}=144.3;
inter.coupling.scalar{H23a,23+numH}=137.2;
inter.coupling.scalar{H16,16+numH}= 146.2;
inter.coupling.scalar{H8,8+numH}=   145.4;
inter.coupling.scalar{H20b,20+numH}=141.0;
inter.coupling.scalar{H20a,20+numH}=141.0;
inter.coupling.scalar{H18a,18+numH}=136.8;
inter.coupling.scalar{H18b,18+numH}=0;
inter.coupling.scalar{H13,13+numH}= 124.4;
inter.coupling.scalar{H17a,17+numH}=133.4;
inter.coupling.scalar{H17b,17+numH}=0;
inter.coupling.scalar{H11a,11+numH}=126.3;
inter.coupling.scalar{H11b,11+numH}=135.9;
inter.coupling.scalar{H14,14+numH}= 130.1;
inter.coupling.scalar{H15b,15+numH}=131.4;
inter.coupling.scalar{H15a,15+numH}=131.4;

% Prune the arrays 
mask=ismember(sys.isotopes,spins);
sys.isotopes=sys.isotopes(mask);
inter.zeeman.scalar=inter.zeeman.scalar(mask);
inter.coordinates=inter.coordinates(mask);
inter.coupling.scalar=inter.coupling.scalar(mask,mask);

end

% Consistency enforcement
function grumble(spins)
if (~iscell(spins))||(~all(cellfun(@ischar,spins)))
    error('spins argument must be a cell array of character strings.');
end
end

% Patience, n. - a minor form of despair, disguised as a virtue.
% 
% Ambrose Bierce

