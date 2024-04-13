% Database of multiplicities and magnetogyric ratios for sta-
% ble and long-lived particles, including spin zero. Syntax:
%
%               [gamma,multiplicity]=spin(name)
%
% Parameters:
%
%      name   - the name of the isotope, e.g. '15N' or
%               '195Pt'; special cases:
%
%                 'G'  - ghost spin: gamma=0, mult=1
%
%                 'N'  - neutron
%
%                 'M'  - muon
%
%                 'P'  - proton
%
%                 'E#' - high-spin electron, # is 
%                        an integer specifying the
%                        multiplicity
%
%                 'C#' - cavity mode, # is an in-
%                        teger specifying the num-
%                        ber of energy levels to
%                        include
%
%                 'V#' - phonon mode, # is an in-
%                        teger specifying the num-
%                        ber of energy levels to
%                        include
%
% Outputs:
%
%      gamma         - magnetogyric ratio, rad/(s*Tesla)
%
%      multiplicity  - spin multiplicity of the particle
%
% Note: data with no source specified was sourced from Google
%       and should be double-checked before running producti-
%       on calculations.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% a.biternas@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spin.m>

function [gamma,multiplicity]=spin(name)

% Check consistency
grumble(name);

% Process the spec
if strcmp(name(1),'E')&&(~isempty(regexp(name,'^E\d','once')))
    
    % High-spin electrons are a special case
    multiplicity=str2double(name(2:end));
    gamma=-1.76085963023e11; % CODATA 2018

elseif strcmp(name(1),'C')&&(~isempty(regexp(name,'^C\d','once')))

    % Cavity modes are a special case
    multiplicity=str2double(name(2:end)); gamma=1;

elseif strcmp(name(1),'V')&&(~isempty(regexp(name,'^V\d','once')))

    % Phonon modes are a special case
    multiplicity=str2double(name(2:end)); gamma=1;
    
else
    
    % Other cases come from the database
    switch name
        case 'G'  % Ghost spin
            multiplicity=1;
            gamma=0;                 % Spin zero particle
        case 'E'  % Electron
            multiplicity=2;
            gamma=-1.76085963023e11; % CODATA 2018
        case 'N'  % Neutron
            multiplicity=2;
            gamma=-1.83247171e8;     % CODATA 2018
        case 'M'  % Muon
            multiplicity=2;
            gamma=-8.51615503e8;     % CODATA 2010
        case {'P','1H'}  % Proton
            multiplicity=2;
            gamma= 2.6752218744e8;   % CODATA 2018
        case '2H'
            multiplicity=3;
            gamma= 4.10662791e7;   % NMR Enc. 1996
        case '3H'
            multiplicity=2;
            gamma= 2.85349779e8;   % NMR Enc. 1996
        case '3He'
            multiplicity=2;
            gamma=-2.03801587e8;   % NMR Enc. 1996
        case '4He'
            multiplicity=1;
            gamma=0;               % Spin 0 nucleus
        case '6Li'
            multiplicity=3;
            gamma= 3.9371709e7;    % NMR Enc. 1996
        case '7Li'
            multiplicity=4;
            gamma= 1.03977013e8;   % NMR Enc. 1996
        case '9Be'
            multiplicity=4;
            gamma=-3.759666e7;     % NMR Enc. 1996
        case '10B'
            multiplicity=7;
            gamma= 2.8746786e7;    % NMR Enc. 1996
        case '11B'
            multiplicity=4;
            gamma= 8.5847044e7;    % NMR Enc. 1996
        case '12C'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '13C'
            multiplicity=2;
            gamma= 6.728284e7;     % NMR Enc. 1996
        case '14N'
            multiplicity=3;
            gamma= 1.9337792e7;    % NMR Enc. 1996
        case '15N'
            multiplicity=2;
            gamma=-2.71261804e7;   % NMR Enc. 1996
        case '16O'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '17O'
            multiplicity=6;
            gamma=-3.62808e7;      % NMR Enc. 1996
        case '18O'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '19F'
            multiplicity=2;
            gamma= 2.518148e8;     % NMR Enc. 1996
        case '20Ne'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '21Ne'
            multiplicity=4;
            gamma=-2.11308e7;      % NMR Enc. 1996
        case '22Ne'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '23Na'
            multiplicity=4;
            gamma=7.0808493e7;     % NMR Enc. 1996
        case '24Mg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '25Mg'
            multiplicity=6;
            gamma=-1.63887e7;      % NMR Enc. 1996
        case '26Mg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '27Al'
            multiplicity=6;
            gamma=6.9762715e7;     % NMR Enc. 1996
        case '28Si'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '29Si'
            multiplicity=2;
            gamma=-5.3190e7;       % NMR Enc. 1996
        case '30Si'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '31P'
            multiplicity=2;
            gamma=10.8394e7;       % NMR Enc. 1996
        case '32S'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '33S'
            multiplicity=4;
            gamma=2.055685e7;      % NMR Enc. 1996
        case '34S'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '36S'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '35Cl'
            multiplicity=4;
            gamma=2.624198e7;      % NMR Enc. 1996
        case '37Cl'
            multiplicity=4;
            gamma=2.184368e7;      % NMR Enc. 1996
        case '36Ar'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '38Ar'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '39Ar'
            multiplicity=8;
            gamma=-1.78e7;         % Needs checking
        case '40Ar'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '39K'
            multiplicity=4;
            gamma= 1.2500608e7;    % NMR Enc. 1996
        case '40K'
            multiplicity=9;
            gamma=-1.5542854e7;    % NMR Enc. 1996
        case '41K'
            multiplicity=4;
            gamma= 0.68606808e7;   % NMR Enc. 1996
        case '40Ca'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '41Ca'
            multiplicity=8;
            gamma=-2.182305e7;     % http://dx.doi.org/10.1103/PhysRevC.63.037301
        case '42Ca'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '43Ca'
            multiplicity=8;
            gamma=-1.803069e7;     % NMR Enc. 1996
        case '44Ca'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '45Sc'
            multiplicity=8;
            gamma=6.5087973e7;     % NMR Enc. 1996
        case '46Ti'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '47Ti'
            multiplicity=6;
            gamma=-1.5105e7;       % NMR Enc. 1996
        case '48Ti'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '49Ti'
            multiplicity=8;
            gamma=-1.51095e7;      % NMR Enc. 1996
        case '50Ti'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '50V'
            multiplicity=13;
            gamma=2.6706490e7;     % NMR Enc. 1996
        case '51V'
            multiplicity=8;
            gamma=7.0455117e7;     % NMR Enc. 1996
        case '52Cr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '53Cr'
            multiplicity=4;
            gamma=-1.5152e7;       % NMR Enc. 1996
        case '54Cr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '55Mn'
            multiplicity=6;
            gamma=6.6452546e7;     % NMR Enc. 1996
        case '54Fe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '56Fe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '57Fe'
            multiplicity=2;
            gamma=0.8680624e7;     % NMR Enc. 1996
        case '58Fe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '59Co'
            multiplicity=8;
            gamma=6.332e7;         % NMR Enc. 1996
        case '60Ni'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '61Ni'
            multiplicity=4;
            gamma=-2.3948e7;       % NMR Enc. 1996
        case '62Ni'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '64Ni'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '63Cu'
            multiplicity=4;
            gamma=7.1117890e7;     % NMR Enc. 1996
        case '65Cu'
            multiplicity=4;
            gamma=7.60435e7;       % NMR Enc. 1996
        case '64Zn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '66Zn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '67Zn'
            multiplicity=6;
            gamma=1.676688e7;      % NMR Enc. 1996
        case '68Zn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '69Ga'
            multiplicity=4;
            gamma=6.438855e7;      % NMR Enc. 1996
        case '71Ga'
            multiplicity=4;
            gamma=8.181171e7;      % NMR Enc. 1996
        case '70Ge'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '72Ge'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '73Ge'
            multiplicity=10;
            gamma=-0.9722881e7;    % NMR Enc. 1996
        case '74Ge'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '75As'
            multiplicity=4;
            gamma=4.596163e7;      % NMR Enc. 1996
        case '74Se'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '76Se'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '77Se'
            multiplicity=2;
            gamma=5.1253857e7;     % NMR Enc. 1996
        case '78Se'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '79Br'
            multiplicity=4;
            gamma=6.725616e7;      % NMR Enc. 1996
        case '81Br'
            multiplicity=4;
            gamma=7.249776e7;      % NMR Enc. 1996
        case '80Kr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '82Kr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '83Kr'
            multiplicity=10;
            gamma=-1.03310e7;      % NMR Enc. 1996
        case '84Kr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '86Kr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '85Rb'
            multiplicity=6;
            gamma=2.5927050e7;     % NMR Enc. 1996
        case '87Rb'
            multiplicity=4;
            gamma=8.786400e7;      % NMR Enc. 1996
        case '84Sr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '86Sr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '87Sr'
            multiplicity=10;
            gamma=-1.1639376e7;    % NMR Enc. 1996
        case '88Sr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '89Y'
            multiplicity=2;
            gamma=-1.3162791e7;    % NMR Enc. 1996
        case '90Zr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '91Zr'
            multiplicity=6;
            gamma=-2.49743e7;      % NMR Enc. 1996
        case '92Zr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '94Zr'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '93Nb'
            multiplicity=10;
            gamma=6.5674e7;        % NMR Enc. 1996
        case '92Mo'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '94Mo'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '95Mo'
            multiplicity=6;
            gamma=-1.7514e7;       % Needs checking
        case '96Mo'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '97Mo'
            multiplicity=6;
            gamma=-1.7884e7;       % Needs checking
        case '98Mo'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '99Tc'
            multiplicity=10;
            gamma=6.0503e7;        % Needs checking
        case '96Ru'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '98Ru'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '99Ru'
            multiplicity=6;
            gamma=-1.2286e7;       % Needs checking
        case '100Ru'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '101Ru'
            multiplicity=6;
            gamma=-1.3773e7;       % Needs checking
        case '102Ru'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '104Ru'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '103Rh'
            multiplicity=2;
            gamma=-0.8468e7;       % NMR Enc. 1996
        case '102Pd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '104Pd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '105Pd'
            multiplicity=6;
            gamma=-1.2305e7;       % Needs checking
        case '106Pd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '108Pd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '110Pd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '107Ag'
            multiplicity=2;
            gamma=-1.0889181e7;    % NMR Enc. 1996
        case '109Ag'
            multiplicity=2;
            gamma=-1.2518634e7;    % NMR Enc. 1996
        case '106Cd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '108Cd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '110Cd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '111Cd'
            multiplicity=2;
            gamma=-5.6983131e7;    % NMR Enc. 1996
        case '112Cd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '113Cd'
            multiplicity=2;
            gamma=-5.9609153e7;    % NMR Enc. 1996
        case '114Cd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '116Cd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '113In'
            multiplicity=10;
            gamma=5.8845e7;        % NMR Enc. 1996
        case '115In'
            multiplicity=10;
            gamma=5.8972e7;        % NMR Enc. 1996
        case '112Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '114Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '115Sn'
            multiplicity=2;
            gamma=-8.8013e7;       % NMR Enc. 1996
        case '116Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '117Sn'
            multiplicity=2;
            gamma=-9.58879e7;      % NMR Enc. 1996
        case '118Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '119Sn'
            multiplicity=2;
            gamma=-10.0317e7;      % NMR Enc. 1996
        case '120Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '122Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '124Sn'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '121Sb'
            multiplicity=6;
            gamma=6.4435e7;        % NMR Enc. 1996
        case '123Sb'
            multiplicity=8;
            gamma=3.4892e7;        % NMR Enc. 1996
        case '120Te'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '122Te'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '123Te' 
            multiplicity=2;
            gamma=-7.059098e7;     % NMR Enc. 1996
        case '124Te'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '125Te'
            multiplicity=2;
            gamma=-8.5108404e7;    % NMR Enc. 1996
        case '126Te'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '128Te'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '130Te'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '127I'
            multiplicity=6;
            gamma=5.389573e7;      % NMR Enc. 1996
        case '124Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '126Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '128Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '129Xe'
            multiplicity=2;
            gamma=-7.452103e7;     % NMR Enc. 1996
        case '130Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '131Xe'
            multiplicity=4;
            gamma=2.209076e7;      % NMR Enc. 1996
        case '132Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '134Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '136Xe'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '133Cs'
            multiplicity=8;
            gamma=3.5332539e7;     % NMR Enc. 1996
        case '130Ba'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '132Ba'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '134Ba'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '135Ba'
            multiplicity=4;
            gamma=2.65750e7;       % NMR Enc. 1996
        case '136Ba'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '137Ba'
            multiplicity=4;
            gamma=2.99295e7;       % NMR Enc. 1996
        case '138Ba'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '138La'
            multiplicity=11;
            gamma=3.557239e7;      % NMR Enc. 1996
        case '139La'
            multiplicity=8;
            gamma=3.8083318e7;     % NMR Enc. 1996
        case '136Ce'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '138Ce'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '139Ce'
            multiplicity=6;
            gamma= 2.906e7;        % Needs checking
        case '140Ce'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '142Ce'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '141Pr'
            multiplicity=6;
            gamma=8.1907e7;        % NMR Enc. 1996
        case '142Nd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '143Nd'
            multiplicity=8;
            gamma=-1.457e7;        % NMR Enc. 1996
        case '144Nd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '145Nd'
            multiplicity=8;
            gamma=-0.898e7;        % NMR Enc. 1996
        case '146Nd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '148Nd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '150Nd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '147Pm'
            multiplicity=8;
            gamma= 3.613e7;        % NMR Enc. 1996
        case '144Sm'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '147Sm'
            multiplicity=8;
            gamma=-1.115e7;        % NMR Enc. 1996
        case '148Sm'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '149Sm'
            multiplicity=8;
            gamma=-0.9192e7;       % NMR Enc. 1996
        case '150Sm'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '152Sm'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '154Sm'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '151Eu'
            multiplicity=6;
            gamma=6.6510e7;        % NMR Enc. 1996
        case '153Eu'
            multiplicity=6;
            gamma=2.9369e7;        % NMR Enc. 1996
        case '152Gd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '154Gd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '155Gd'
            multiplicity=4;
            gamma=-0.82132e7;      % NMR Enc. 1996
        case '156Gd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '157Gd'
            multiplicity=4;
            gamma=-1.0769e7;       % NMR Enc. 1996
        case '158Gd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '160Gd'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '159Tb'
            multiplicity=4;
            gamma=6.4306e7;        % Needs checking
        case '156Dy'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '158Dy'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '160Dy'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '161Dy'
            multiplicity=6;
            gamma=-0.9201e7;       % NMR Enc. 1996
        case '162Dy'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '163Dy'
            multiplicity=6;
            gamma=1.289e7;         % NMR Enc. 1996
        case '164Dy'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '165Ho'
            multiplicity=8;
            gamma=5.710e7;         % NMR Enc. 1996
        case '162Er'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '164Er'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '166Er'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '167Er'
            multiplicity=8;
            gamma=-0.77157e7;      % NMR Enc. 1996
        case '168Er'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '170Er'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '169Tm'
            multiplicity=2;
            gamma=-2.218e7;        % NMR Enc. 1996
        case '168Yb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '170Yb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '171Yb'
            multiplicity=2;
            gamma=4.7288e7;        % NMR Enc. 1996
        case '172Yb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '173Yb'
            multiplicity=6;
            gamma=-1.3025e7;       % NMR Enc. 1996
        case '174Yb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '176Yb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '175Lu'
            multiplicity=8;
            gamma=3.0552e7;        % NMR Enc. 1996
        case '176Lu'
            multiplicity=15;
            gamma=2.1684e7;        % NMR Enc. 1996
        case '174Hf'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '176Hf'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '177Hf'
            multiplicity=8;
            gamma=1.086e7;         % NMR Enc. 1996
        case '178Hf'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '179Hf'
            multiplicity=10;
            gamma=-0.6821e7;       % NMR Enc. 1996
        case '180Hf'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '180Ta'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '181Ta'
            multiplicity=8;
            gamma=3.2438e7;        % NMR Enc. 1996
        case '180W'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '182W'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '183W'
            multiplicity=2;
            gamma=1.1282403e7;     % NMR Enc. 1996
        case '184W'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '186W'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '185Re'
            multiplicity=6;
            gamma=6.1057e7;        % NMR Enc. 1996
        case '187Re'
            multiplicity=6;
            gamma=6.1682e7;        % NMR Enc. 1996
        case '184Os'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '186Os'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '187Os'
            multiplicity=2;
            gamma=0.6192895e7;     % NMR Enc. 1996
        case '188Os'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '189Os'
            multiplicity=4;
            gamma=2.10713e7;       % NMR Enc. 1996
        case '190Os'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '192Os'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '191Ir'
            multiplicity=4;
            gamma=0.4812e7;        % NMR Enc. 1996
        case '193Ir'
            multiplicity=4;
            gamma=0.5227e7;        % NMR Enc. 1996
        case '190Pt'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '192Pt'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '194Pt'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '195Pt'
            multiplicity=2;
            gamma=5.8385e7;        % NMR Enc. 1996
        case '196Pt'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '198Pt'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '197Au'
            multiplicity=4;
            gamma=0.473060e7;      % NMR Enc. 1996
        case '196Hg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '198Hg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '199Hg'
            multiplicity=2;
            gamma=4.8457916e7;     % NMR Enc. 1996
        case '200Hg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '201Hg'
            multiplicity=4;
            gamma=-1.788769e7;     % NMR Enc. 1996
        case '202Hg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '204Hg'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '203Tl'
            multiplicity=2;
            gamma=15.5393338e7;    % NMR Enc. 1996
        case '205Tl'
            multiplicity=2;
            gamma=15.6921808e7;    % NMR Enc. 1996
        case '204Pb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '206Pb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '207Pb'
            multiplicity=2;
            gamma=5.58046e7;       % NMR Enc. 1996
        case '208Pb'
            multiplicity=1;        
            gamma=0;               % Spin 0 nucleus
        case '209Bi'
            multiplicity=10;
            gamma=4.3750e7;        % NMR Enc. 1996
        case '209Po'
            multiplicity=2;
            gamma= 7.4e7;          % NMR Enc. 1996
        case '227Ac'
            multiplicity=4;
            gamma= 3.5e7;          % NMR Enc. 1996
        case '229Th'
            multiplicity=6;
            gamma= 0.40e7;         % NMR Enc. 1996
        case '231Pa'
            multiplicity=4;
            gamma=3.21e7;          % NMR Enc. 1996
        case '235U'
            multiplicity=8;
            gamma=-0.4926e7;       % Needs checking
        case '237Np'
            multiplicity=6;
            gamma=3.1e7;           % NMR Enc. 1996
        case '239Pu'
            multiplicity=2; 
            gamma=0.972e7;         % NMR Enc. 1996
        case '241Am' 
            multiplicity=6;
            gamma=1.54e7;          % NMR Enc. 1996
        case '243Am' 
            multiplicity=6;
            gamma=1.54e7;          % NMR Enc. 1996
        case '247Cm'
            multiplicity=10;
            gamma=0.20e7;          % NMR Enc. 1996
        case {'210At', '222Rn', '223Fr', '226Ra', '247Bk', '251Cf', '252Es', '253Es', '257Fm', '258Md', '259No'}
            error([name ' - no data available in the current NMR literature.']);
        otherwise
            error([name ' - unknown isotope.']);
    end
    
end

end

% Consistency enforcement
function grumble(name)
if ~ischar(name)
    error('isotope specification must be a character string.');
end
end

% I mean that there is no way to disarm any man except through guilt.
% Through that which he himself has accepted as guilt. If a man has 
% ever stolen a dime, you can impose on him the punishment intended
% for a bank robber and he will take it. He'll bear any form of mise-
% ry, he'll feel that he deserves no better. If there's not enough 
% guilt in the world, we must create it. If we teach a man that it's
% evil to look at spring flowers and he believes us and then does it,
% we'll be able to do whatever we please with him. He won't defend 
% himself. He won't feel he's worth it. He won't fight. But save us
% from the man who lives up to his own standards. Save us from the 
% man of clean conscience. He's the man who'll beat us.
%
% Ayn Rand, "Atlas Shrugged"

