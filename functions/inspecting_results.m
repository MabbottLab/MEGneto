function inspecting_results(paths)

load([paths.anout_grp '/fcp_5_powspctrm.mat']);
atlas = importdata('/mnt/sda/juanita/fieldtrip/template/atlas/aal/ROI_MNI_V4.txt');
disp('Press any key to show next AAL region.')

for region = 1:90
    imagesc(squeeze(mean(pow_spctrm(:,region,:,:),1)))
    xlabel('Timewindows')
    ylabel('Frequencies 1-100Hz')
    title(sprintf('%s power spectrum', atlas.textdata{region,2}),'Interpreter','none');
    pause;
end