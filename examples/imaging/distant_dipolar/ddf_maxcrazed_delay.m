% Unfitted MAXCRAZED delay prediction for Chen et al., Fig. 5
% (doi:10.1063/1.4721637). The three finite lobes after the mixing
% pulse and the two free intervals occupy 60 ms in total. The first
% free interval separates the strong lobes; the second follows them.
% Outputs: inter-lobe gaps (s), four normalised peak heights and
% three hetero/homo ratios per gap, and both complex FID dimensions.

function result=ddf_maxcrazed_delay()

% Reserve the time occupied by the 2X and two strong finite lobes
free_time=60e-3-3*1.5e-3;
result.delay=[0 15 30 45 55.5]*1e-3;
for k=1:numel(result.delay)

    % Move the return lobe without changing time from mixing to readout
    data=ddf_maxcrazed_scan([12 16 16],0.25e-3,result.delay(k),...
                            max(0,free_time-result.delay(k)));
    result.fid(:,:,:,k)=data.fid;
    result.peaks(:,:,k)=data.predicted;
    result.ratios(k,:)=sum(data.predicted(2:3,:),1)./sum(data.predicted([1 4],:),1);

end

end


