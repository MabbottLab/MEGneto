function [res] = detectionBadChannels_Algorithm( data, thr )
       
res1 = zeros(size(data.label,1),size(data.trial,2));
res2 = zeros(size(data.label,1),size(data.trial,2));
checkData = data.trial{1,1};
sizeTrialData   = size(checkData,1);
        
% standardize each channel
for indexe = 1:sizeTrialData          
    channelData    = checkData(indexe,:);
    channelData(isnan(channelData)) = [];
    channelMedian    = median(channelData,2,'omitnan'); %median
    channelStd     = std(channelData,0,2,'omitnan'); %std 
    channelStdAll(indexe) = (channelStd/sqrt(sizeTrialData)); %standard error
    channelData(isnan(channelData)) = [];
    channelNormalized   = ( channelData - channelMedian ) / channelStd;    
    trialNormalized(indexe,:) = channelNormalized;
end;

trialMedian = median(trialNormalized,1,'omitnan');
chanNorms = 1:sizeTrialData;
meanSTD = mean(channelStdAll);
    
if ~isempty(trialMedian) % if the last segment happens to be all NaNs

    for indexe = 1:sizeTrialData
        channelDataNorm     = trialNormalized(indexe,:);  
        channelNorm         = norm(trialMedian - channelDataNorm);  % no NaNs allowed 
        chanNorms(1,indexe) = channelNorm; 
        channelNormalizedSTD = abs( meanSTD - channelStdAll(indexe))/meanSTD;
        trialSTDNormized(indexe) = channelNormalizedSTD;
    end;    


    % figure;    plot(chanNorms);
    % figure; plot(channelStdAll);
    % figure; plot(trialSTDNormized);
    [F,XI] = ksdensity(chanNorms);
    % [F,XI] = ksdensity(trialSTDNormized); % Don't use
    % figure;  plot(XI,F);   
    [~,maxIndex] = max(F);
    [Fs,XIs] = ksdensity(channelStdAll); 
    % figure;  plot(Fs);
    [~,maxIndexs] = max(Fs);


    if maxIndex <= 50 % less than 50% 
        pthr = prctile(XI,thr);
        [~,indT] = min(abs(XI-pthr));     
        cutOff = XI(1,indT);        
%         cutOff = XI(1,maxIndex*2); 
    else
         cutOff = XI(1,thr);
    end;
    disp(cutOff);
    for indexe = 1:sizeTrialData
        if chanNorms(1,indexe) > cutOff 
            res1(indexe,1) = 1;                        
        end;
    end;

    if maxIndexs <= 50
        cutOff2 = XIs(1,maxIndexs*2); % less than 0.1%
    else
        cutOff2 = XIs(1,thr);
    end
    disp(cutOff2);
    for indexe = 1:sizeTrialData
        if channelStdAll(indexe) > cutOff2 
            res2(indexe,1) = 1;                        
        end;

    end;
end    
    
    res = res1.*res2; % needs to be found in mean and std search
    
end

        
