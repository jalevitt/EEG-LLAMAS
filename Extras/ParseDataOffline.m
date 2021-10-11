
Paths = {'C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19.mat',...
    'C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19_jeff.mat',...
    'C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19_jeff2.mat'...
    };

WindowsPerSecond = 2;
TimePoints = {[291:1/WindowsPerSecond:320, 351:1/WindowsPerSecond:380], [321:1/WindowsPerSecond:350, 381:1/WindowsPerSecond:410];
              [41:1/WindowsPerSecond:100, 161:1/WindowsPerSecond:220], [101:1/WindowsPerSecond:160, 221:1/WindowsPerSecond:280];
              [111:1/WindowsPerSecond:170, 231:1/WindowsPerSecond:290], [171:1/WindowsPerSecond:230, 291:1/WindowsPerSecond:350]
}; 
EC = [];
EO = [];

Channel = 22;
for F=1:length(Paths)
    EEG = LoadRealTimeFile(Paths{F});
    data = EEG.data(Channel, :);
    data = data(1, 1:find(data, 1, 'last'));
    l = length(data);
    fs = EEG.fs;
    alphaVec = zeros(1000, 1);
    thetaVec = zeros(1000, 1);
    betaVec = zeros(1000, 1);
    smoothedAlphaVec = zeros(1000, 1);
    smoothedThetaVec = zeros(1000, 1);
    smoothedBetaVec = zeros(1000, 1);
    
    highPassIIR = designfilt('highpassiir', 'FilterOrder', 4, 'PassbandFrequency', 5, 'SampleRate', EEG.fs);
    [psd, bins] = pwelch(rand(EEG.fs,1), [], [], [], EEG.fs);
    [~, fullMin] = min(abs(bins - 3));
    [~, fullMax] = min(abs(bins - 30));
    [~, alphaMin] = min(abs(bins - 8));
    [~, alphaMax] = min(abs(bins - 12));
    [~, thetaMin] = min(abs(bins - 4));
    [~, thetaMax] = min(abs(bins - 7));
    [~, betaMin] = min(abs(bins - 13));
    [~, betaMax] = min(abs(bins - 30));
    sawTooth = linspace(0,1,5);
    %sawTooth = [0, 0, 1, 1, 1];
    %sawTooth = hamming(5);
    time = fs * 5 + 1;
    
    while time <= l - fs
        
        sample = data(time:time + fs/WindowsPerSecond - 1);
        FiltData = filtfilt(highPassIIR, sample);
        
        [psd, ~] = pwelch(FiltData, [], [], [], EEG.fs);

        psd = psd/sum(psd(fullMin:fullMax));
        
        idx = (time - 1)/(EEG.fs/WindowsPerSecond);
        alphaVec(idx) = mean(psd(alphaMin:alphaMax));
        smoothedAlphaVec(idx) = dot(alphaVec(idx - 4:idx), sawTooth);
        
        thetaVec(idx) = mean(psd(thetaMin:thetaMax));
        smoothedThetaVec(idx) = dot(thetaVec(idx - 4:idx), sawTooth);

        betaVec(idx) = mean(psd(betaMin:betaMax));
        smoothedBetaVec(idx) = dot(betaVec(idx - 4:idx), sawTooth);

        
        time = time + fs/WindowsPerSecond;
    end
    EC = [EC; [smoothedAlphaVec(TimePoints{F,1} * WindowsPerSecond), smoothedThetaVec(TimePoints{F,1} * WindowsPerSecond), smoothedBetaVec(TimePoints{F,1} * WindowsPerSecond)]];
    EO = [EO; [smoothedAlphaVec(TimePoints{F,2} * WindowsPerSecond), smoothedThetaVec(TimePoints{F,2} * WindowsPerSecond), smoothedBetaVec(TimePoints{F,2} * WindowsPerSecond)]];
end
%%
figure
histogram(EC(:, 1), 0:0.01:.5)
hold on
histogram(EO(:, 1), 0:0.01:.5)
hold off

Labels = [ones(length(EC), 1); zeros(length(EO), 1)];
Features = [EC(:, 1); EO(:, 1)];
%Features = [Features, ones(size(Features))];
t = templateNaiveBayes();
%t = templateSVM('Standardize', 'on');
%EOECModel = fitcnb(Features, Labels, 'CrossVal', 'on', 'KFold', 10);
EOECModel = fitcecoc(Features, Labels, 'KFold', 10, 'Learners', t);
[Predictions, score] = kfoldPredict(EOECModel);
Accuracy = 100 * sum(Predictions == Labels) / length(Labels);
Accuracy
[X, Y, T, AUC] = perfcurve(Labels, score(:, 2), 1);
AUC
figure
plot(X, Y)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Receveiver Operating Characteristic')

EOECModel = fitcecoc(Features, Labels, 'Learners', t);