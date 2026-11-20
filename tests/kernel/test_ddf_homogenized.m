% Check the HOMOGENIZED spectrum against spin echoes integrated over a
% real-space solvent kernel, and check the acquired line positions. Syntax:
%
%                    result=test_ddf_homogenized()
%
% Outputs:
%
%    result - Spinach regression test result structure

function result=test_ddf_homogenized()

% Run the actual example and restore the path on return
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path=path;
restore_path=onCleanup(@()path(old_path));
addpath(fullfile(root_dir,'examples','imaging','distant_dipolar'));
data=ddf_homogenized();
result=new_test_result('kernel/ddf_homogenized','Solvent-referenced spectroscopy',...
                       'Resolve a J doublet and measure finite-distance broadening.');

% Read the two strongest resolved indirect lines from each echo trace
expect=[116 124;116 124;116 124;116 124;132 140];
for c=1:5

    % Sum complex direct bins to recover the first point of direct acquisition
    trace=abs(sum(data.spectrum(:,:,c),2));
    peaks=find(islocalmax(trace));
    [~,idx]=maxk(trace(peaks),2);
    found=sort(-data.f1(peaks(idx)));
    result=test_close(result,sprintf('component positions, case %d',c),...
                      found,expect(c,:)',0.5,0,'within one native indirect bin');

end

% Common offsets must shift the direct axis, not merely disappear
result=test_close(result,'direct frequency shift',data.spectrum(:,:,4),...
                  circshift(data.spectrum(:,:,2),8,2),1e-13,0,...
                  '16 Hz common offset moves the direct spectrum by eight bins');
result=test_close(result,'reversed ramp',data.spectrum(:,:,2),...
                  data.spectrum(:,:,3),1e-13,0,'symmetric illuminated volume');

% Broad conventional spectrum bounded by the same single-pulse calibration
spread=sum(abs(data.sq).^2,2).^2./sum(abs(data.sq).^4,2);
result=test_true(result,'single quantum broadening',spread(2)>4*spread(1),...
                 'more than four times as many effectively occupied frequency bins');
result=test_true(result,'linewidth growth',all(diff(data.width)>0)&&...
                 data.width(end)>4*data.width(1),'isolated echo component FWHM');

% Prepare one commensurate 3D spatial period for an independent voxel-field echo
npts=[16 16 16];
[x,y,z]=ndgrid((0:15)*2*pi/16);
grating=cos(x+y+z);
magn=[zeros(numel(x),2) grating(:)];

% Numerically average explicit hard-pulse rotations at distinct grating phases
fid=ifft2(ifftshift(data.spectrum(:,:,2)))*numel(data.spectrum(:,:,2));
azimuth=2*pi*(0:127)/128;
beta=pi/4;
k=pi/data.distance(1);
for idx=[1 39 141]

    % Integrate the exponential real-space field kernel in half-pitch units
    time=data.t1(idx);
    field=integral(@(u)exp(-u).*cos(2*pi*2e4*time*u/k),0,Inf)/3;
    signal=0;
    for freq=[116 124]

        % Rotate the excited solute through mixing and conjugate at the refocuser
        phase=azimuth+2*pi*freq*time;
        transverse=sin(phase)+1i*cos(beta)*cos(phase);
        precession=exp(-1i*0.02*field*sin(beta)*cos(azimuth));
        signal=signal+mean(transverse.*precession)/2;

    end

    % Compare the oblique solvent mode with the actual Cartesian field backend
    if time>0

        % Make each Cartesian component commensurate without rotating the dipolar axis
        wave=[k 2*pi*2e4*time*sqrt(2/3) 2*pi*2e4*time/sqrt(3)];
        spacing=2*pi./(npts.*wave);
        symbol=ddf_periodic(npts,spacing,1/3);
        voxel_field=ddf_cartesian(magn,npts,symbol);
        result=test_close(result,sprintf('3D field, t1 index %d',idx),...
                          voxel_field(:,3)/1.25663706127e-6,field*grating(:),...
                          1e-10,0,'Cartesian voxel field versus real-space integral');

    end

    % Compare complex echo amplitudes without fitting receiver phase or gain
    result=test_close(result,sprintf('real-space pulse oracle, t1 index %d',idx),...
                      fid(idx,1),exp(-pi*time)*signal,1e-12,0,...
                      'real-space field integral and explicit rotations');

end

end


