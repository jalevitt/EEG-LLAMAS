
Paths = {'C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19.mat',...
    'C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19_jeff.mat',...
    'C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19_jeff2.mat'...
    };

fs = 200
Channel = 22;
WindowsPerSecond = 4;

TimePoints = {[291:1/WindowsPerSecond:320, 351:1/WindowsPerSecond:380], [321:1/WindowsPerSecond:350, 381:1/WindowsPerSecond:410];
              [41:1/WindowsPerSecond:100, 161:1/WindowsPerSecond:220], [101:1/WindowsPerSecond:160, 221:1/WindowsPerSecond:280];
              [111:1/WindowsPerSecond:170, 231:1/WindowsPerSecond:290], [171:1/WindowsPerSecond:230, 291:1/WindowsPerSecond:350]
}; 

PredictionWindowSize = fs;
Features = zeros(1, fs / WindowsPerSecond);
Labels = zeros(1, PredictionWindowSize);

RawSignals = zeros(1, fs / WindowsPerSecond);
FiltSignals = zeros(1, fs / WindowsPerSecond);
RawTargets = zeros(1, PredictionWindowSize);
Mansouris = zeros(1, PredictionWindowSize);
PhaseBins = zeros(1, 19);
PLV  = [];
UseEyesOpen = 0;
sampleNumber = 1;
for F=1:length(Paths)
    EEG = LoadRealTimeFile(Paths{F});
    data = EEG.data(Channel, :);
    data = data(1, 1:find(data, 1, 'last'));
    l = length(data);
    %data = 3 * sin(linspace(0, (l - 1) / fs, l) * 11.4 * 2 * pi) + randn(1, l);

    fs = EEG.fs;
    %bandPassIIR = designfilt('bandpassfir', 'FilterOrder', 4, 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'SampleRate', EEG.fs);
    bandPassIIR = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 15,  'SampleRate', EEG.fs);
    bandPassIIR2 = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 50,  'SampleRate', EEG.fs);
    
    time = fs * 5 + 1;
    temp = filtfilt(bandPassIIR2, data);
    phase = angle(hilbert(temp - mean(temp)));
    while time < l - 2 * fs
        if sum((time - 1)/fs == TimePoints{F, 1}) > 0 || UseEyesOpen
            sample = data(time:time + fs/WindowsPerSecond - 1);
            RawSignals(sampleNumber, :) = (sample - mean(sample)) / std(sample);
            target = data(time + fs/WindowsPerSecond:time + fs/WindowsPerSecond + PredictionWindowSize - 1);
            RawTargets(sampleNumber, :) = (target - mean(target)) / std(target);

            FiltData = filtfilt(bandPassIIR, sample);
            FiltSignals(sampleNumber, :) = FiltData;
            FiltTarget = filtfilt(bandPassIIR, target);
            M = MansouriPhase(FiltData, fs, PredictionWindowSize);
            compLag = 40 + round(rand * 5);
            [~, trigger] = min(M(1 + compLag:end)); 
            triggerIDX = time + fs/WindowsPerSecond - 1 + trigger + compLag;
            triggerPhase = phase(triggerIDX);
            p = round(9.499 * triggerPhase/pi) + 10;
            PhaseBins(p) = PhaseBins(p) + 1;
            
            PLV(sampleNumber) = PhaseLockValue(M(1:10), RawTargets(sampleNumber, 1:10));
            Mansouris(sampleNumber, :) = M;

            sampleNumber = sampleNumber + 1;
        end
        time = time + fs/WindowsPerSecond;
    end
end

figure
bar(PhaseBins)

%%
figure
histogram(PLV)
xlabel('Phase Lock Value')
ylabel('Count')
title('True Phase Lock Histogram')

P = randperm(sampleNumber - 1);
randMansouris = Mansouris(P, :);
randPLV = zeros(sampleNumber - 1, 1);
for i = 1:(sampleNumber - 1)
    randPLV(i) = PhaseLockValue(RawTargets(i, 1:10), randMansouris(i, 1:10));
end

figure
histogram(randPLV)
xlabel('Phase Lock Value')
ylabel('Count')
title('Randomized Phase Lock Histogram')

%%
figure
A = 0;
TrueVec = 0;
subplot(2, 1, 1)
hold on
for n = 1:(sampleNumber - 1)

    Ha = angle(hilbert(RawTargets(n, 1:10)));
    Hb = angle(hilbert(Mansouris(n, 1:10)));

    RoseA = exp(1).^(1i * (Ha - Hb));

    R = mean(RoseA);
    TrueVec = TrueVec + R;
  
    line([0, real(R)], [0, imag(R)])
    xlim([-1, 1])
    ylim([-1, 1])

end
hold off
A = A / (sampleNumber - 1);

TrueVec = TrueVec / (sampleNumber - 1);

subplot(2, 1, 2)
hold on
RandVec = 0;
for n = 1:(sampleNumber - 1)

    Ha = angle(hilbert(RawTargets(n, 1:10)));
    Hc = angle(hilbert(randMansouris(n, 1:10)));

    RoseB = exp(1).^(1i * (Ha - Hc));
    
    R = mean(RoseB);
    RandVec = RandVec + R;
  
    line([0, real(R)], [0, imag(R)])
    xlim([-1, 1])
    ylim([-1, 1])
end
RandVec = RandVec / (sampleNumber -1);
hold off

%%

PLV2 = zeros(sampleNumber - 1, PredictionWindowSize);
Theta = zeros(sampleNumber - 1, PredictionWindowSize);
randPLV2 = zeros(sampleNumber - 1, PredictionWindowSize);
randTheta = zeros(sampleNumber - 1, PredictionWindowSize);
for n = 1:PredictionWindowSize
    for i = 1:sampleNumber - 1
        [randPLV2(i, n), randTheta(i, n) ]= PhaseLockValue(RawTargets(i, 1:n), randMansouris(i, 1:n));
        [PLV2(i, n), Theta(i, n)] = PhaseLockValue(RawTargets(i, 1:n), Mansouris(i, 1:n));
    end
end

M = mean(PLV2);
MeanAngleError = mean(abs(Theta));
randM = mean(randPLV2);
randMeanAngleError = mean(abs(randTheta));

figure
subplot(2,2,1)
plot(1:PredictionWindowSize, M)
xlabel('Prediction Length')
ylabel('Phase Lock Value')
title('True Values')
ylim([0, 1])

subplot(2, 2, 2)
plot(1:PredictionWindowSize, MeanAngleError)
hold on
line([0, PredictionWindowSize], [pi/2, pi/2], 'Color', 'k')
hold off
xlabel('Prediction Length')
ylabel('Mean Absolute Phase Error (Radians)')
title('True Values')
ylim([0, pi])

subplot(2, 2, 3)
plot(1:PredictionWindowSize, randM)
xlabel('Prediction Length')
ylabel('Phase Lock Value')
title('Randomized Values')
ylim([0, 1])

subplot(2, 2, 4)
plot(1:PredictionWindowSize, randMeanAngleError)
hold on
line([0, PredictionWindowSize], [pi/2, pi/2], 'Color', 'k')
hold off
xlabel('Prediction Length')
ylabel('Mean Absolute Phase Error (Radians)')
title('Randomized Values')
ylim([0, pi])