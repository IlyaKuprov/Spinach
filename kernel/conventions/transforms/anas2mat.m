% Converts anisotropy and asymmetry representation of a 3x3 interaction
% tensor (Haeberlen-Mehring convention) into the corresponding matrix.
% Euler angles should be specified in radians. Syntax:
%
%                 M=anas2mat(iso,an,as,alp,bet,gam)
%
% Parameters:
%
%        iso  - isotropic part of the interaction, defined as
%               (xx+yy+zz)/3 in terms of eigenvaues
%
%         an  - interaction anisotropy, defined as zz-(xx+yy)/2
%               in terms of eigenvalues
%
%         as  - interaction asymmetry, defined as (yy-xx)/(zz-iso)
%               in terms of eigenvalues
%
%        alp  - alpha Euler angle in radians
%
%        bet  - beta Euler angle in radians
%
%        gam  - gamma Euler angle in radians
%
% Outputs:
%
%          M  - 3x3 matrix
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=anas2mat.m>

function M=anas2mat(iso,an,as,alp,bet,gam)

% Check consistency
grumble(iso,an,as,alp,bet,gam);

% Compute reduced anisotropy
ra=2*an/3;

% Compute eigenvalues
zz=iso+ra;
yy=iso-ra*(1-as)/2;
xx=iso-ra*(1+as)/2;

% Rotate the matrix
R=euler2dcm(alp,bet,gam);
M=R*diag([xx yy zz])*R';

end

% Consistency enforcement
function grumble(iso,an,as,alp,bet,gam)
if (~isnumeric(iso))||(~isreal(iso))||(~isscalar(iso))||...
   (~isnumeric(an))||(~isreal(an))||(~isscalar(an))||...
   (~isnumeric(as))||(~isreal(as))||(~isscalar(as))||...
   (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('all inputs must be real scalars.');
end
end

% Why is it that beauty is no longer the standard by which we judge things,
% asked a female friend of mine as I was writing this column. [...] I told
% her why. Beauty is elitist, carried off by a few, and we now live in the
% age of the common and vulgar. Beauty is an optimistic feeling, and we ha-
% ve all now become pessimists. Beauty means great inequality, and that's 
% the real no-no of our times. We are supposed to all be equal now, at least
% that's what the crooks who rule us and make these ludicrous laws try to
% force us to believe.
%
% But we are not. The busybodies think they can cure the ills of an unequal
% society by keeping the smart ones back. That's an old trick that was tri-
% ed with comprehensive schools and one that ended in utter failure. [...]
% No wonder opinion-makers loathe beauty, it is as unequal as it gets. Be-
% auty does not need to promote itself over others, it's too self-evident.
%
% And I forgot about ugly art. The sanctification of Picasso, not to menti-
% on the grotesque Lucian Freud, are perfect examples of visual ugliness 
% triumphing over the sublime art of Leonardo, Rembrandt, Sargent and Degas.
% [...] Critics now judge art not on aesthetics but on whether it advances
% a progressive political agenda. When I see the crap that sells for hund-
% reds of millions I really want to shout.
%
% Taki Theodoracopulos

