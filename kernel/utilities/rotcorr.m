% Rigid-molecule rotational diffusion from a solvent-aware contact-surface
% ellipsoid. Syntax:
%
%   [tau,D,axes_len]=rotcorr(atom_symbols,xyz,solvent,temperature)
%
% Parameters:
%
%   atom_symbols - N-by-1 cell array of element character row vectors;
%                  case-sensitive supported symbols: H B C N O F Si P S Cl Se Br I
%   xyz          - real finite floating-point N-by-3 Cartesian atom
%                  coordinates, Angstrom; no coincident atom centres
%   solvent      - character row vector: 'water' or 'chloroform'
%   temperature  - absolute temperature, Kelvin; 273.16-373 K (water),
%                  210-334 K (chloroform)
%
% Outputs:
%
%   tau       - isotropic-equivalent rank-2 time, 1/(2*trace(D)), s
%   D         - 3-by-3 rotational diffusion tensor in the input
%               Cartesian frame, inverse seconds
%   axes_len  - 3-by-1 ellipsoid semi-axes, descending, Angstrom
%
% Mantina elemental van der Waals radii are used in both solvents.
% Water adds the ROTDIF3 2.2 Angstrom shell and 1.4 Angstrom probe.
% Chloroform uses an unhydrated envelope (zero shell and zero probe).
% Pure H2O and CHCl3 viscosities are temperature dependent at about
% ambient pressure; D2O and CDCl3 are not represented.
% Parameters and their sources are shipped in etc/rotcorr_sources.md.
% Use complete atom lists for one rigid body in dilute pure solvent.
% No isotope substitution, mixtures, pressure correction, flexibility,
% aggregation, or specific solvent binding is modelled.
%
% Exposed contact points on atomic spheres are area weighted; surface
% covariance gives semi-axes sqrt(3*eigenvalues), and Perrin integrals
% give diffusion. This ELM-inspired approximation is not SURF: re-entrant
% patches are omitted. Protein hydration is not a small-molecule fit;
% ellipsoidal reduction is unreliable for strongly non-globular bodies.
% Chloroform predictions are uncalibrated continuum estimates, not a
% validated solvation model. Numerical convergence is not model accuracy.
%
% Fibonacci directions per atom double from 256 until tau, all semi-axes,
% and diffusion in every direction change by less than 1% on two successive
% refinements. This is a relative self-consistency test, not an absolute
% error bound. Failure at 1048576 directions per atom raises an error;
% this cap bounds the number of rows in quadrature work arrays.
% Each refinement costs O(N^2+N*npoints*neighbours), with O(N+npoints)
% workspace. Strongly occluded or irregular bodies may be expensive.
%
% Based on Ryabov et al., JACS 128 (2006), 15432-15444,
% https://doi.org/10.1021/ja062715t. The scalar tau is not a universal
% single-exponential time for anisotropic rotation; use D for that case.
%
% talos@spindynamics.org
%
% <https://spindynamics.org/wiki/index.php?title=rotcorr.m>

function [tau,D,axes_len]=rotcorr(atom_symbols,xyz,solvent,temperature)

% Check consistency
grumble(atom_symbols,xyz,solvent,temperature);

% Select the documented atomic and solvent model
symbols={'H';'B';'C';'N';'O';'F';'Si';'P';'S';'Cl';'Se';'Br';'I'};
radii=[1.10;1.92;1.70;1.55;1.52;1.47;2.10;1.80;1.80;1.75;1.90;1.83;1.98];
[known,indices]=ismember(atom_symbols,symbols);
if ~all(known)
    error('unsupported element symbol; see the rotcorr header.');
end
shell=radii(indices); temperature=double(temperature);
switch solvent
    case 'water'
        shell=shell+2.2; probe=1.4; reduced_temp=temperature/300;
        visc=1e-6*(280.68*reduced_temp^(-1.9)+511.45*reduced_temp^(-7.7)+...
                  61.131*reduced_temp^(-19.6)+0.45903*reduced_temp^(-40));
    case 'chloroform'
        probe=0;
        visc=exp(-14.109+1049.2/temperature+0.5377*log(temperature));
end

% Remove an irrelevant translation before accumulating moments
xyz=double(xyz); xyz=xyz-mean(xyz,1);
if any(~isfinite(xyz),'all')
    error('centred coordinates exceed floating-point range.');
end

% Double the angular resolution until all outputs stabilise twice
previous_tau=Inf; previous_D=nan(3); previous_axes=inf(3,1);
stable=0;
for npoints=2.^(8:20)
    [tau,D,axes_len]=surface_estimate(xyz,shell,probe,...
                                    temperature,visc,npoints);
    if ~isfinite(tau)||any(~isfinite(D),'all')||any(eig(D)<=0)
        error('geometry produces nonfinite or nonpositive diffusion.');
    end
    if isfinite(previous_tau)
        changes=[abs(tau-previous_tau)/tau;...
                 abs(axes_len-previous_axes)./axes_len;...
                 abs(eig(D-previous_D,D))];
        if max(changes)<0.01
            stable=stable+1;
        else
            stable=0;
        end
        if stable==2, return; end
    end
    previous_tau=tau; previous_D=D; previous_axes=axes_len;
end
error('surface sampling did not converge to 1%% at 1048576 points per atom.');

end

% Evaluate the contact-surface ellipsoid at one angular resolution
function [tau,D,axes_len]=surface_estimate(xyz,shell,probe,temperature,visc,npoints)

% Build an equal-solid-angle Fibonacci sphere
z=1-2*((1:npoints)'-0.5)/npoints;
phi=pi*(3-sqrt(5))*(0:npoints-1)';
dirs=[sqrt(1-z.^2).*cos(phi) sqrt(1-z.^2).*sin(phi) z];

% Accumulate area-weighted exposed surface moments
area=0; first=zeros(1,3); second=zeros(3);
for n=1:size(xyz,1)

    % Only overlapping expanded spheres can occlude this atom
    delta=xyz-xyz(n,:); dist=sqrt(sum(delta.^2,2));
    neighbours=find((dist<shell+shell(n)+2*probe)&...
                    ((1:size(xyz,1))'~=n));
    centres=xyz(n,:)+(shell(n)+probe)*dirs;
    exposed=true(npoints,1);
    for k=neighbours'
        delta=centres-xyz(k,:);
        exposed=exposed&(sum(delta.^2,2)>=(shell(k)+probe)^2);
    end

    % Integrate contact-surface moments with sphere-area weights
    points=xyz(n,:)+shell(n)*dirs(exposed,:);
    weight=shell(n)^2/npoints;
    area=area+weight*size(points,1);
    first=first+weight*sum(points,1);
    second=second+weight*(points'*points);
end

% Obtain the surface-covariance ellipsoid
if area==0
    error('no exposed surface samples at the current resolution.');
end
centre=first/area; covar=second/area-centre'*centre;
if any(~isfinite(covar),'all')
    error('surface covariance exceeds floating-point range.');
end
[frame,values]=eig(covar,'vector');
[values,order]=sort(values,'descend'); frame=frame(:,order);
if any(~isfinite(values))||any(values<=0)
    error('surface covariance is nonfinite or singular.');
end
axes_len=sqrt(3*values);

% Evaluate dimensionless Perrin integrals to avoid SI underflow
scale=max(axes_len); squares=(axes_len/scale).^2;
ints=zeros(3,1);
for n=1:3
    ints(n)=integral(@(u)1./((squares(n)+u).*...
        sqrt((squares(1)+u).*(squares(2)+u).*(squares(3)+u))),...
        0,Inf,'RelTol',1e-10,'AbsTol',0);
end

% Convert ellipsoid friction to diffusion in the original frame
terms=squares.*ints;
rates=3*1.380649e-23*temperature/(16*pi*visc*(scale*1e-10)^3)*...
      (sum(terms)-terms)./(sum(squares)-squares);
D=frame*diag(rates)*frame'; tau=1/(2*sum(rates));

end

% Validate the chemical, geometric, and solvent specifications
function grumble(atom_symbols,xyz,solvent,temperature)
if ~isfloat(xyz)||~isreal(xyz)||size(xyz,2)~=3||...
   ~ismatrix(xyz)||isempty(xyz)||any(~isfinite(xyz),'all')
    error('xyz must be a nonempty real finite N-by-3 array.');
end
if ~iscell(atom_symbols)||~isequal(size(atom_symbols),[size(xyz,1) 1])||...
   ~all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),atom_symbols))
    error('atom_symbols must be an N-by-1 cell array of character row vectors.');
end
if size(unique(xyz,'rows'),1)~=size(xyz,1)
    error('atom coordinates must be distinct.');
end
if ~ischar(solvent)||~isrow(solvent)||...
   ~ismember(solvent,{'water','chloroform'})
    error('solvent must be ''water'' or ''chloroform''.');
end
if ~isnumeric(temperature)||~isreal(temperature)||~isscalar(temperature)||...
   ~isfinite(temperature)
    error('temperature must be a real finite numeric scalar in Kelvin.');
end
if strcmp(solvent,'water')&&(temperature<273.16||temperature>373)
    error('water temperature must be between 273.16 and 373 Kelvin.');
end
if strcmp(solvent,'chloroform')&&(temperature<210||temperature>334)
    error('chloroform temperature must be between 210 and 334 Kelvin.');
end
end


