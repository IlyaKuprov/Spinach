% Spin system of allyl pyruvate. Isotropic chemical shifts and 
% J-couplings determined by spectral fitting, coordinates and
% chemical shift anisotropies from DFT. Syntax:
%
%                [sys,inter]=allyl_pyruvate(spins)
%
% Parameters:
%
%    spins   - a cell array containing the isotopes
%              to import, e.g. {'1H','13C'}
%
% Outputs:
%
%    sys, inter - Spinach data structures with the
%                 specification of the spin system
%
% Note: 13C-13C J-couplings are not provided - this spin system
%       is for natural abundance 13C simulations only.
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=allyl_pyruvate.m>

function [sys,inter]=allyl_pyruvate(spins)

% Check consistency
grumble(spins);

% Spin system
sys.isotopes={'13C','13C','13C','13C','13C','13C',...
              '1H','1H','1H','1H','1H','1H','1H','1H'};
sys.labels={'C1','C2','C3','C4','C5','C6',...
            'Ha','Hb','Hc','Hd1','Hd2','He1','He2','He3'};

% Chemical shifts (from fitting)
inter.zeeman.scalar={119.9643 130.6962 66.8463 160.4202 191.6305 ...
                      26.6961   5.4086  5.3273   5.9678   4.7466 ...
                       4.7466   2.4836  2.4836   2.4836};

% Chemical shift anisotropies (DFT)
inter.zeeman.matrix{idxof(sys,'C1')}=-remtrace([ 73.8088   56.3125   -5.7848
                                                 71.4477   44.7750   82.9998
                                                  2.6990   82.7441   67.5866]);
inter.zeeman.matrix{idxof(sys,'C2')}=-remtrace([ 65.5131   69.8026    0.8618
                                                 58.8550   40.5199   92.9966
                                                 -0.5172   90.4267   56.9090]);
inter.zeeman.matrix{idxof(sys,'C3')}=-remtrace([140.4211   17.8082    1.7258
                                                  9.3520  110.9204   -5.9851
                                                  3.3228   -1.5838  104.2861]);
inter.zeeman.matrix{idxof(sys,'C4')}=-remtrace([-64.9706  -30.5508   -9.0450
                                                -15.1976   74.5683   -0.9381
                                                 -9.1799   -1.2291   81.5741]);
inter.zeeman.matrix{idxof(sys,'C5')}=-remtrace([-83.4746  -39.6016  -10.1866
                                                -47.0307  -33.5938   -2.3508
                                                 -9.9814   -2.2044   95.1185]);
inter.zeeman.matrix{idxof(sys,'C6')}=-remtrace([184.8168   -3.1537    1.8807
                                                -13.0905  139.2936   -0.7330
                                                  1.9002   -0.1601  152.4313]);
inter.zeeman.matrix{idxof(sys,'Ha')}=-remtrace([ 28.2575    1.3410   -3.0000
                                                 -0.2827   23.6667    2.7125
                                                 -1.8783    1.3010   25.2580]);
inter.zeeman.matrix{idxof(sys,'Hb')}=-remtrace([ 28.0342   -0.1603   -1.2955
                                                  2.1711   23.8330    0.8358
                                                 -2.6492    1.2849   25.8721]);
inter.zeeman.matrix{idxof(sys,'Hc')}=-remtrace([ 27.9364   -0.1637   -1.6581
                                                 -1.2265   23.5562    1.4542
                                                  0.3094    0.3715   24.3166]);
inter.zeeman.matrix{idxof(sys,'Hd1')}=-remtrace([28.0196    1.6136    1.4936
                                                 -0.1493   26.5655    2.4262
                                                  0.9592    3.0435   26.1123]);
inter.zeeman.matrix{idxof(sys,'Hd2')}=-remtrace([29.6956    1.5093   -1.5351
                                                  1.0802   26.5728   -2.4805
                                                 -2.1842   -2.0104   25.5956]);
inter.zeeman.matrix{idxof(sys,'He1')}=-remtrace([29.3456   -0.8213   -1.3863
                                                  1.1157   28.3559    1.4973
                                                 -1.7101    3.3511   28.9577]);
inter.zeeman.matrix{idxof(sys,'He2')}=-remtrace([28.9925   -0.4100    1.3818
                                                  1.3247   28.3222   -1.4380
                                                  1.6857   -3.4936   29.3033]);
inter.zeeman.matrix{idxof(sys,'He3')}=-remtrace([32.9096    4.0247    0.3478
                                                  0.1631   29.6606   -0.0268
                                                  0.3721    0.2046   26.7796]);

% 13C-1H J-couplings (from fitting)
inter.coupling.scalar=cell(14,14);
inter.coupling.scalar{idxof(sys,'C1'),idxof(sys,'Ha')}= 156.04; 
inter.coupling.scalar{idxof(sys,'C1'),idxof(sys,'Hb')}= 160.11; 
inter.coupling.scalar{idxof(sys,'C1'),idxof(sys,'Hd1')}=  5.57; 
inter.coupling.scalar{idxof(sys,'C1'),idxof(sys,'Hd2')}=  5.57; 
inter.coupling.scalar{idxof(sys,'C2'),idxof(sys,'Hc')}= 158.96;
inter.coupling.scalar{idxof(sys,'C2'),idxof(sys,'Hd1')}= -4.44;
inter.coupling.scalar{idxof(sys,'C2'),idxof(sys,'Hd2')}= -4.44;
inter.coupling.scalar{idxof(sys,'C2'),idxof(sys,'Ha')}=  -3.16; 
inter.coupling.scalar{idxof(sys,'C2'),idxof(sys,'Hb')}=  -0.85; 
inter.coupling.scalar{idxof(sys,'C3'),idxof(sys,'Hd1')}=149.00;
inter.coupling.scalar{idxof(sys,'C3'),idxof(sys,'Hd2')}=149.00;
inter.coupling.scalar{idxof(sys,'C3'),idxof(sys,'Hc')}= -13.89;
inter.coupling.scalar{idxof(sys,'C3'),idxof(sys,'Ha')}=   4.62; 
inter.coupling.scalar{idxof(sys,'C3'),idxof(sys,'Hb')}=   8.04; 
inter.coupling.scalar{idxof(sys,'C4'),idxof(sys,'Hd1')}=  3.21;
inter.coupling.scalar{idxof(sys,'C4'),idxof(sys,'Hd2')}=  3.21;
inter.coupling.scalar{idxof(sys,'C4'),idxof(sys,'He1')}=  1.60;
inter.coupling.scalar{idxof(sys,'C4'),idxof(sys,'He2')}=  1.60;
inter.coupling.scalar{idxof(sys,'C4'),idxof(sys,'He3')}=  1.60;
inter.coupling.scalar{idxof(sys,'C5'),idxof(sys,'He1')}= -6.68;
inter.coupling.scalar{idxof(sys,'C5'),idxof(sys,'He2')}= -6.68;
inter.coupling.scalar{idxof(sys,'C5'),idxof(sys,'He3')}= -6.68;
inter.coupling.scalar{idxof(sys,'C6'),idxof(sys,'He1')}=129.52;
inter.coupling.scalar{idxof(sys,'C6'),idxof(sys,'He2')}=129.52;
inter.coupling.scalar{idxof(sys,'C6'),idxof(sys,'He3')}=129.52;

% 1H-1H J-couplings (from fitting)
inter.coupling.scalar{idxof(sys,'Ha'),idxof(sys,'Hb')}=   +1.16;
inter.coupling.scalar{idxof(sys,'Ha'),idxof(sys,'Hc')}=  +17.18;
inter.coupling.scalar{idxof(sys,'Hb'),idxof(sys,'Hc')}=  +10.40;
inter.coupling.scalar{idxof(sys,'Ha'),idxof(sys,'Hd1')}=  -1.52;
inter.coupling.scalar{idxof(sys,'Ha'),idxof(sys,'Hd2')}=  -1.52;
inter.coupling.scalar{idxof(sys,'Hb'),idxof(sys,'Hd1')}=  -1.16;
inter.coupling.scalar{idxof(sys,'Hb'),idxof(sys,'Hd2')}=  -1.16;
inter.coupling.scalar{idxof(sys,'Hc'),idxof(sys,'Hd1')}=  +5.98;
inter.coupling.scalar{idxof(sys,'Hc'),idxof(sys,'Hd2')}=  +5.98;

% Coordinates (from DFT)
inter.coordinates={[-3.845334  -0.310277  -0.426714]   % C1
                   [-2.765389  -0.342326   0.350942]   % C2
                   [-1.608016   0.574585   0.183489]   % C3
                   [ 0.723669   0.439738   0.024494]   % C4
                   [ 1.939407  -0.502779  -0.043635]   % C5
                   [ 3.252146   0.196785  -0.121022]   % C6
                   [-4.696729  -0.972873  -0.251664]   % Ha
                   [-3.921746   0.391176  -1.264816]   % Hb
                   [-2.698234  -1.050632   1.185688]   % Hc
                   [-1.494575   1.257208   1.043829]   % Hd1
                   [-1.697969   1.193186  -0.724943]   % Hd2
                   [ 3.374780   0.856946   0.752099]   % He1
                   [ 3.273104   0.855016  -1.003844]   % He2
                   [ 4.068811  -0.533089  -0.167358]}; % He3

% Decide which atoms to return
mask=ismember(sys.isotopes,spins);

% Prune the arrays 
sys.isotopes=sys.isotopes(mask);
sys.labels=sys.labels(mask);
inter.zeeman.scalar=inter.zeeman.scalar(mask);
inter.zeeman.matrix=inter.zeeman.matrix(mask);
inter.coordinates=inter.coordinates(mask);
inter.coupling.scalar=inter.coupling.scalar(mask,mask);

end

% Consistency enforcement
function grumble(spins)
if (~iscell(spins))||(~all(cellfun(@ischar,spins)))
    error('spins argument must be a cell array of character strings.');
end
end

% You will sag in your chains as embers scorch your
% feet. It will be endless. The word 'endless' - do
% you understand it?
%
% Vladimir Bulgakov

