function [mvldata] = mvl_calc(LFsigtemp,HFsigtemp)
% calculate  mean vector length (complex value) per trial
% mvldata dim: LF*HF
    LFphas   = angle(LFsigtemp);
    HFamp    = abs(HFsigtemp);
    mvldata = nanmean(HFamp(1,:).*exp(1i*LFphas(1,:)));
end