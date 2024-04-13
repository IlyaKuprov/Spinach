% Spin system of cyprinol. Isotropic chemical shifts and J-couplings
% are taken from http://dx.doi.org/10.1002/mrc.4782 and, when not gi-
% ven there, estimated by tossing a twenty-sided coin. Syntax:
%
%                    [sys,inter,bas]=cyprinol()
%
% Outputs:
%
%    sys   - Spinach spin system description structure
%
%    inter - Spinach interaction description structure
%
%    bas   - Spinach basis set description structure
%
% Note: if you are looking for a test spin system, strychnine.m is a
%       more complete alternative.
%
% Bud MacAulay
% Ilya Kuprov
%
% <https://spindynamics.org/wiki/index.php?title=cyprinol.m>

function [sys,inter,bas]=cyprinol()

% Hydrogen atoms
H_iso=cell(1,42); H_iso(:)={'1H'}; numH=numel(H_iso);

% Shorthands for human-readable coupling designations below
H1a=1;   H1b=2;   H2=2;    H3=3;    H4a=5;   H4b=6;   H5=7;    H6a=8;   
H6b=9;   H7=10;   H8=11;   H9=12;   H11a=13; H11b=14; H12=15;  H14=16;  
H15a=17; H15b=18; H16a=19; H16b=20; H17=21;  H18a=22; H19a=23; H20=24;  
H21a=25; H22a=26; H22b=27; H23a=28; H23b=29; H24a=30; H24b=31; H25=32;  
H26a=33; H26b=34; H27a=35; H27b=36; H18b=37; H18c=38; H19b=39; H19c=40; 
H21b=41; H21c=42;

% Proton chemical shifts
H_CS=cell(1,numH);
H_CS{H1a}= 1.42;  H_CS{H17}= 1.85;  H_CS{H1b}= 1.42;  H_CS{H18a}=0.72;
H_CS{H2}=  1.64;  H_CS{H19a}=0.82;  H_CS{H3}=  3.99;  H_CS{H20}= 1.40;
H_CS{H4a}= 1.51;  H_CS{H21a}=1.02;  H_CS{H4b}= 1.32;  H_CS{H22a}=1.46;
H_CS{H5}=  2.15;  H_CS{H22b}=1.11;  H_CS{H6a}= 1.42;  H_CS{H23a}=1.48;
H_CS{H6b}= 1.35;  H_CS{H23b}=1.31;  H_CS{H7}=  3.79;  H_CS{H24a}=1.35;
H_CS{H8}=  1.48;  H_CS{H24b}=1.35;  H_CS{H9}=  1.68;  H_CS{H25}= 1.80;
H_CS{H11a}=1.68;  H_CS{H26a}=4.05;  H_CS{H11b}=1.57;  H_CS{H26b}=4.00;
H_CS{H12}= 3.95;  H_CS{H27a}=3.60;  H_CS{H14}= 1.94;  H_CS{H27b}=3.56;
H_CS{H15a}=1.77;  H_CS{H18b}=0.72;  H_CS{H15b}=1.12;  H_CS{H18c}=0.72;
H_CS{H16a}=1.89;  H_CS{H19b}=0.82;  H_CS{H16b}=1.28;  H_CS{H19c}=0.82;
H_CS{H21b}=1.02;  H_CS{H21c}=1.02;

% TODO: proton J-couplings

% TODO: proton coordinates

% Carbon atoms
C_iso=cell(1,27); C_iso(:)={'13C'}; numC=numel(C_iso);

% Carbon chemical shifts
C_CS=cell(1,numC);
C_CS{1}= 31.71;  C_CS{15}=22.72;  C_CS{2}= 28.09;  C_CS{16}=27.31;
C_CS{3}= 65.73;  C_CS{17}=47.08;  C_CS{4}= 35.13;  C_CS{18}=11.45;
C_CS{5}= 31.24;  C_CS{19}= 8.99;  C_CS{6}= 36.25;  C_CS{20}=35.76;
C_CS{7}= 67.35;  C_CS{21}=16.74;  C_CS{8}= 39.81;  C_CS{22}=35.98;
C_CS{9}= 39.02;  C_CS{23}=23.25;  C_CS{10}=35.50;  C_CS{24}=27.84;
C_CS{11}=27.80;  C_CS{25}=40.91;  C_CS{12}=72.61;  C_CS{26}=67.56;
C_CS{13}=46.04;  C_CS{27}=61.73;  C_CS{14}=41.73;

% Carbon J-coupling estimates
CC_J=cell(numC);
CC_J{1,2}=40;    CC_J{1,10}=40;   CC_J{1,3}=3;
CC_J{1,5}=3;     CC_J{1,4}=0.3;   CC_J{1,6}=0.3;
CC_J{2,3}=40;    CC_J{2,4}=3;     CC_J{2,10}=3;
CC_J{2,5}=0.3;   CC_J{2,9}=0.3;   CC_J{3,4}=40;
CC_J{3,5}=3;     CC_J{3,10}=0.3;  CC_J{3,6}=0.3;
CC_J{4,5}=40;    CC_J{4,6}=3;     CC_J{4,10}=3;
CC_J{4,7}=0.3;   CC_J{4,9}=0.3;   CC_J{5,6}=40;
CC_J{5,10}=40;   CC_J{5,7}=3;     CC_J{5,9}=3;
CC_J{5,11}=0.3;  CC_J{5,8}=0.3;   CC_J{6,7}=40;
CC_J{6,8}=3;     CC_J{6,10}=3;    CC_J{6,9}=0.3;
CC_J{7,8}=40;    CC_J{7,9}=3;     CC_J{7,14}=3;
CC_J{7,10}=0.3;  CC_J{7,11}=0.3;  CC_J{7,13}=0.3;
CC_J{7,15}=0.3;  CC_J{8,9}=40;    CC_J{8,14}=40;
CC_J{8,10}=3;    CC_J{8,11}=3;    CC_J{8,15}=3;
CC_J{8,13}=3;    CC_J{8,16}=0.3;  CC_J{8,12}=0.3;
CC_J{8,17}=0.3;  CC_J{9,10}=40;   CC_J{9,11}=40;
CC_J{9,12}=3;    CC_J{9,14}=3;    CC_J{9,13}=0.3;
CC_J{9,15}=0.3;  CC_J{11,12}=40;  CC_J{11,13}=3;
CC_J{11,14}=0.3; CC_J{11,17}=0.3; CC_J{12,13}=40;
CC_J{12,14}=3;   CC_J{12,17}=3;   CC_J{12,15}=0.3;
CC_J{12,16}=0.3; CC_J{12,20}=0.3; CC_J{13,14}=40;
CC_J{13,17}=40;  CC_J{13,15}=3;   CC_J{13,20}=3;
CC_J{14,15}=40;  CC_J{14,16}=3;   CC_J{14,17}=3;
CC_J{14,20}=0.3; CC_J{15,16}=40;  CC_J{15,17}=3;
CC_J{15,20}=0.3; CC_J{16,17}=40;  CC_J{16,20}=3;
CC_J{16,21}=0.3; CC_J{16,22}=0.3; CC_J{17,20}=40;
CC_J{17,21}=3;   CC_J{17,22}=3;   CC_J{17,23}=0.3;
CC_J{20,21}=40;  CC_J{20,22}=40;  CC_J{20,23}=3;
CC_J{20,24}=0.3; CC_J{21,22}=3;   CC_J{21,23}=0.3;
CC_J{22,23}=40;  CC_J{22,24}=3;   CC_J{22,25}=0.3;
CC_J{23,24}=40;  CC_J{23,25}=3;   CC_J{23,26}=0.3;
CC_J{23,27}=0.3; CC_J{24,25}=40;  CC_J{24,26}=3;
CC_J{24,27}=3;   CC_J{25,26}=40;  CC_J{25,27}=40;
CC_J{26,27}=3;   CC_J{18,13}=40;  CC_J{18,14}=3;
CC_J{18,12}=3;   CC_J{18,17}=3;   CC_J{18,11}=0.3;
CC_J{18,16}=0.3; CC_J{18,8}=0.3;  CC_J{18,15}=0.3;
CC_J{19,10}=40;  CC_J{19,1}=3;    CC_J{19,2}=0.3;
CC_J{19,9}=3;    CC_J{19,5}=3;    CC_J{19,4}=0.3;
CC_J{19,6}=0.3;  CC_J{19,8}=0.3;
 
% TODO: carbon coordinates

% Carbon-proton J-coupling estimates
CH_J=cell(numC,numH);
CH_J{1,H1a}=150;   CH_J{1,H1b}=150;   CH_J{2,H2}=150;    CH_J{4,H4a}=150;
CH_J{4,H4b}=150;   CH_J{5,H5}=150;    CH_J{6,H6a}=150;   CH_J{6,H6b}=150;
CH_J{7,H7}=150;    CH_J{8,H8}=150;    CH_J{9,H9}=150;    CH_J{11,H11a}=150;
CH_J{11,H11b}=150; CH_J{12,H12}=150;  CH_J{14,H14}=150;  CH_J{15,H15a}=150;
CH_J{15,H15b}=150; CH_J{16,H16a}=150; CH_J{16,H16b}=150; CH_J{17,H17}=150;
CH_J{18,H18a}=150; CH_J{18,H18b}=150; CH_J{18,H18c}=150; CH_J{19,H19a}=150;
CH_J{19,H19b}=150; CH_J{19,H19c}=150; CH_J{20,H20}=150;  CH_J{21,H21a}=150;
CH_J{21,H21b}=150; CH_J{21,H21c}=150; CH_J{22,H22a}=150; CH_J{22,H22b}=150;
CH_J{23,H23a}=150; CH_J{23,H23b}=150; CH_J{24,H24a}=150; CH_J{24,H24b}=150;
CH_J{25,H25}=150;  CH_J{26,H26a}=150; CH_J{26,H26b}=150; CH_J{27,H27a}=150;
CH_J{27,H27b}=150;

CH_J{1,H2}=40;     CH_J{2,H1a}=40;    CH_J{2,H1b}=40;    CH_J{2,H3}=40;
CH_J{3,H2}=40;     CH_J{3,H4a}=40;    CH_J{3,H4b}=40;    CH_J{4,H3}=40;
CH_J{4,H5}=40;     CH_J{5,H4a}=40;    CH_J{5,H4b}=40;    CH_J{5,H6a}=40;
CH_J{5,H6b}=40;    CH_J{6,H5}=40;     CH_J{6,H7}=40;     CH_J{7,H6a}=40;
CH_J{7,H6b}=40;    CH_J{7,H8}=40;     CH_J{8,H7}=40;     CH_J{8,H9}=40;
CH_J{8,H14}=40;    CH_J{9,H8}=40;     CH_J{9,H11a}=40;   CH_J{9,H11b}=40;
CH_J{10,H1a}=40;   CH_J{10,H1b}=40;   CH_J{10,H19a}=40;  CH_J{10,H19b}=40;
CH_J{10,H19c}=40;  CH_J{11,H9}=40;    CH_J{11,H12}=40;
CH_J{12,H11a}=40;  CH_J{12,H11b}=40;  CH_J{13,H14}=40;   CH_J{13,H17}=40;
CH_J{13,H18a}=40;  CH_J{13,H18b}=40;  CH_J{13,H18c}=40;  CH_J{14,H8}=40;
CH_J{14,H15a}=40;  CH_J{14,H15b}=40;  CH_J{15,H14}=40;   CH_J{15,H16a}=40;
CH_J{15,H16b}=40;  CH_J{16,H15a}=40;  CH_J{16,H15b}=40;  CH_J{16,H17}=40;
CH_J{17,H16a}=40;  CH_J{17,H16b}=40;  CH_J{17,H20}=40;   CH_J{20,H21a}=40;
CH_J{20,H21b}=40;  CH_J{20,H21c}=40;  CH_J{20,H22a}=40;  CH_J{20,H22b}=40;
CH_J{20,H17}=40;   CH_J{21,H20}=40;   CH_J{22,H20}=40;   CH_J{22,H23a}=40;
CH_J{22,H23b}=40;  CH_J{23,H22a}=40;  CH_J{23,H22b}=40;  CH_J{23,H24a}=40;
CH_J{23,H24b}=40;  CH_J{24,H23a}=40;  CH_J{24,H23b}=40;  CH_J{24,H25}=40;
CH_J{25,H24a}=40;  CH_J{25,H24b}=40;  CH_J{25,H26a}=40;  CH_J{25,H26b}=40;
CH_J{25,H27a}=40;  CH_J{25,H27b}=40;  CH_J{26,H25}=40;   CH_J{27,H25}=40;

CH_J{1,H3}=2;      CH_J{1,H5}=2;      CH_J{1,H9}=2;      CH_J{1,H19a}=2;
CH_J{1,H19b}=2;    CH_J{1,H19c}=2;    CH_J{3,H5}=2;      CH_J{3,H1a}=2;
CH_J{3,H1b}=2;     CH_J{5,H1a}=2;     CH_J{5,H1b}=2;     CH_J{5,H3}=2;
CH_J{5,H7}=2;      CH_J{5,H9}=2;      CH_J{7,H5}=2;      CH_J{7,H9}=2;
CH_J{8,H11a}=2;    CH_J{8,H11b}=2;    CH_J{8,H15a}=2;    CH_J{8,H15b}=2;
CH_J{9,H1a}=2;     CH_J{9,H1b}=2;     CH_J{9,H5}=2;      CH_J{9,H7}=2;
CH_J{9,H12}=2;     CH_J{9,H14}=2;     CH_J{9,H19a}=2;    CH_J{9,H19b}=2;
CH_J{9,H19c}=2;    CH_J{10,H6a}=2;    CH_J{10,H6b}=2;    CH_J{10,H11a}=2;
CH_J{10,H11b}=2;   CH_J{11,H8}=2;     CH_J{12,H9}=2;     CH_J{12,H14}=2;
CH_J{12,H18a}=2;   CH_J{12,H18b}=2;   CH_J{12,H18c}=2;   CH_J{13,H11a}=2;
CH_J{13,H11a}=2;   CH_J{13,H16a}=2;   CH_J{13,H16b}=2;   CH_J{14,H9}=2;
CH_J{14,H12}=2;    CH_J{14,H18a}=2;   CH_J{14,H18b}=2;   CH_J{14,H18c}=2;
CH_J{15,H8}=2;     CH_J{17,H14}=2;    CH_J{17,H18a}=2;   CH_J{17,H18b}=2;
CH_J{17,H18c}=2;   CH_J{17,H21a}=2;   CH_J{17,H21b}=2;   CH_J{17,H21c}=2;
CH_J{17,H22a}=2;   CH_J{17,H22b}=2;   CH_J{18,H14}=2;    CH_J{18,H12}=2;
CH_J{18,H17}=2;    CH_J{19,H1a}=2;    CH_J{19,H1b}=2;    CH_J{19,H5}=2;
CH_J{19,H9}=2;     CH_J{20,H23a}=2;   CH_J{20,H23b}=2;   CH_J{21,H17}=2;
CH_J{21,H22a}=2;   CH_J{21,H22b}=2;   CH_J{22,H17}=2;    CH_J{22,H21a}=2;
CH_J{22,H21b}=2;   CH_J{22,H21c}=2;   CH_J{23,H20}=2;    CH_J{23,H25}=2;
CH_J{24,H22a}=2;   CH_J{24,H22b}=2;   CH_J{24,H26a}=2;   CH_J{24,H26b}=2;
CH_J{24,H27a}=2;   CH_J{24,H27b}=2;   CH_J{25,H23a}=2;   CH_J{25,H23b}=2;
CH_J{26,H24a}=2;   CH_J{26,H24b}=2;   CH_J{26,H27a}=2;   CH_J{26,H27b}=2;
CH_J{27,H24a}=2;   CH_J{27,H24b}=2;   CH_J{27,H26a}=2;   CH_J{27,H26b}=2;

% Combine isotope arrays
sys.isotopes=[H_iso C_iso];

% Combine chemical shift arrays
inter.zeeman.scalar=[H_CS C_CS];

% Combine J-coupling arrays
inter.coupling.scalar=cell(numH+numC);
inter.coupling.scalar((numH+1):end,(numH+1):end)=CC_J;
inter.coupling.scalar((numH+1):end,1:numH)=CH_J;

% Symmetry settings
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[H21a H21b H21c],...
               [H18a H18b H18c],...
               [H19a H19b H19c]};

end

% My husband and I are either going to buy a dog or
% have a child. We can't decide whether to ruin our
% carpets or ruin our lives.
% 
% Rita Rudner

