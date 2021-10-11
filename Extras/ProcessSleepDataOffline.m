load('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\Test_11_20_19_Sleep1.mat');
%%
%{
BandPassFIR = designfilt('bandpassfir', ...
    'PassbandFrequency1', 0.5, ...
    'Passbandfrequency2', 2, ...
    'StopbandFrequency1', 0.1 ,...
    'StopbandFrequency2', 10, ...
    'StopbandAttenuation1', 15, ...
    'StopbandAttenuation2', 10, ...
    'PassbandRipple', 1,  ...
    'SampleRate', EEG.fs);
%}

BandPassFIR = designfilt('bandpassiir', ...
    'PassbandFrequency1', 0.5, ...
    'Passbandfrequency2', 2, ...
    'StopbandFrequency1', 0.1 ,...
    'StopbandFrequency2', 10, ...
    'StopbandAttenuation1', 20, ...
    'StopbandAttenuation2', 20, ...
    'PassbandRipple', 0.1, ...
    'DesignMethod', 'butter', ...
    'SampleRate', EEG.fs);
BandPassFIR2 = designfilt('bandpassfir', ...
    'PassbandFrequency1', 0.5, ...
    'Passbandfrequency2', 2, ...
    'StopbandFrequency1', 0.1 ,...
    'StopbandFrequency2', 10, ...
    'StopbandAttenuation1', 60, ...
    'StopbandAttenuation2', 60, ...
    'PassbandRipple', 0.1, ...
    'SampleRate', EEG.fs);
HighPassFIR = designfilt('highpassfir', ...
    'PassbandFrequency', 0.4, ...
    'StopbandFrequency', 0.1 ,...
    'StopbandAttenuation', 40, ...
    'PassbandRipple', 1, ...
    'SampleRate', EEG.fs);
[b, a] = tf(BandPassFIR);

%%
PrimaryChannel = 18 + 32; %Fz

Signal = filter(b, a, EEG.Recording(:, PrimaryChannel));
Signal2 = filtfilt(BandPassFIR2, EEG.Recording(:, PrimaryChannel));

n = length(Signal);
Time = linspace(1/EEG.fs, n / EEG.fs, n);

figure
plot(Time/60, Signal2)
xlabel('Time (s)')
ylabel('voltage')
hold on
plot(Time/60, EEG.Recording(:, PrimaryChannel))
plot(Time / 60, Signal)
legend({'Strong Filter', 'No Filter', 'Weak Filter'})

%% METHOD A: Mansouri Prediction
figure
hold on

WindowsPerSecond = 4;
Threshold = 20;
StimTimes = [];
count = 1;
FilteredSignal = zeros(n, 1);
[temp, z] = filter(b, a, EEG.Recording(1:time - 1, PrimaryChannel));
FilteredSignal(1:time - 1) = temp;

time = 1 + 1000 * EEG.fs;
while time < n - EEG.fs / WindowsPerSecond
    sample = EEG.Recording(time:time + EEG.fs / WindowsPerSecond - 1, PrimaryChannel);
    [filtSample, z] = filter(b, a, sample, z);
    FilteredSignal(time:time+EEG.fs / WindowsPerSecond - 1) = filtSample;
    if max(filtSample) > Threshold
        [Prediction] = MansouriPhase(FilteredSignal(time - EEG.fs * 3:time+EEG.fs / WindowsPerSecond - 1)', EEG.fs, EEG.fs * 2);
        [~, offset] = max(Prediction);
        StimTimes(count) = time + WindowsPerSecond - 1 + offset;
        count = count + 1;
    end
    
    time = time + EEG.fs / WindowsPerSecond;
end


hold on
for i = 1:length(StimTimes)
    line([StimTimes(i)/EEG.fs, StimTimes(i)/EEG.fs], [-40, 40], 'Color', [1, 0, 0]);
end
plot(Time, FilteredSignal)
hold off

Phase = [];
H = hilbert(Signal);
for i = 1:length(StimTimes)
    Phase = [Phase, angle(H(StimTimes(i)))];
end

figure
histogram(Phase, linspace(-1 * pi, pi, 18));
xlabel('Phase (radians)')
ylabel('Hits')

%% Method B: Immediate Response
figure
hold on
%time = 1 + 1000 * EEG.fs;
time = 1;
Threshold = 30.1;
StimTimes = [];
LastStim = 0;
buffer = EEG.fs;
count = 1;
z = zeros(6, 1);
FilteredSignal = zeros(n, 1);
WindowsPerSecond = 4;
%[temp, z] = filter(b, a, EEG.Recording(1:time - 1, PrimaryChannel));
%FilteredSignal(1:time - 1) = temp;

while time < n - EEG.fs / WindowsPerSecond
    StepSize = round(rand * 4) + 1;
    sample = EEG.Recording(time:time + StepSize - 1, PrimaryChannel);
    [filtSample, z] = filter(b, a, sample, z);
    FilteredSignal(time:time+StepSize - 1) = filtSample;
    if max(filtSample) > Threshold && ((time + StepSize - buffer) > LastStim)
        StimTimes(count) = time + StepSize - 1 + 3;
        LastStim = time + StepSize - 1 + 3;
        count = count + 1;
    end
    time = time + StepSize;
end


for i = 1:length(StimTimes)
    line([StimTimes(i)/EEG.fs, StimTimes(i)/EEG.fs], [-40, 40], 'Color', [1, 0, 0]);
end
plot(Time, FilteredSignal)
%plot(Time, Signal2)
hold off

%%
Phase = [];
H = hilbert(Signal2);
for i = 1:length(StimTimes)
    Phase = [Phase, angle(H(StimTimes(i)))];
end

figure
histogram(Phase, linspace(-1 * pi, pi, 18));
xlabel('Phase (radians)')
ylabel('Hits')

%%
[NumericalSleepScore, ~, RawSleepScore] = xlsread('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\rtEEG_score2.csv');
%%
Counts.N1 = 0;
Counts.N2 = 0;
Counts.N3 = 0;
Counts.Wake = 0;
Duration.N1 = 0;
Duration.N2 = 0;
Duration.N3 = 0;
Duration.Wake = 0;

for i = 1:length(RawSleepScore) - 1
    state = RawSleepScore{i, 3};
    d = NumericalSleepScore(i, 2) - NumericalSleepScore(i, 1);
    Duration.(state) = Duration.(state) + d;
    for stim = 1:length(StimTimes)
        if StimTimes(stim)/(EEG.fs * 60) < NumericalSleepScore(i, 2) && StimTimes(stim)/(EEG.fs * 60) > NumericalSleepScore(i, 1)
            Counts.(state) = Counts.(state) + 1;
        end
    end
    
end

figure
bar([Duration.Wake, Duration.N1, Duration.N2, Duration.N3])
xticklabels({'Wake', 'N1', 'N2', 'N3'})
ylabel('Duration (min)')
title('Time in each state')

figure
bar([Counts.Wake / Duration.Wake, Counts.N1 / Duration.N1, Counts.N2 / Duration.N2, Counts.N3 / Duration.N3])
xticklabels({'Wake', 'N1', 'N2', 'N3'})
title('Hit rate in each state')
ylabel('Hits / Minute')

%%
[n, chans] = size(EEG.Recording);
StimTimes = 1:n;
StimTimes(EEG.Recording(:, chans - 1) ~= 32958) = [];

Phase = [];
H = hilbert(Signal2);
for i = 1:length(StimTimes)
    Phase = [Phase, angle(H(StimTimes(i) + 8))];
end

figure
histogram(Phase, linspace(-1 * pi, pi, 18));
xlabel('Phase (radians)')
ylabel('Hits')

%%
n = length(NumericalSleepScore);
t = 0:0.5:NumericalSleepScore(n - 1,2) - 0.5;
PSD = zeros(length(t), 513);
labels = zeros(length(t), 1);
colors = zeros(length(t), 3);
ColorKey.Wake = [1, 0, 0];
ColorKey.N1 = [0, 1, 0];
ColorKey.N2 = [0, 0, 1];
ColorKey.N3 = [0.5, 0.5, 0.5];
LabelKey.Wake = 1;
LabelKey.N1 = 2;
LabelKey.N2 = 3;
LabelKey.N3 = 4;
currenttp = 1;

figure
hold on
for i = 1:length(t)
    startTime = t(i) * 60 * EEG.fs + 1;
    endTime = startTime + 30 * EEG.fs - 1;
    tpSample = EEG.Recording(startTime:endTime, PrimaryChannel);
    filtered_tpSample = filtfilt(HighPassFIR, tpSample);
    [PSDtemp, freqBins] = pwelch(filtered_tpSample, hamming(EEG.fs * 4), EEG.fs/2, [], EEG.fs);
    PSD(i, :) = PSDtemp';
    if t(i) == NumericalSleepScore(currenttp + 1, 1)
        currenttp = currenttp + 1;
    end
    state = RawSleepScore{currenttp, 3};
    colors(i, :) = ColorKey.(state);
    labels(i) = LabelKey.(state);
    loglog(freqBins, PSDtemp, 'Color', ColorKey.(state))
end
hold off
xlabel('frequency (Hz)')
ylabel('Power (uV^2)')

[coeff, score, latent, tsquared, explained] = pca(PSD);
figure
hold on
for i = 1:3
    plot(freqBins, coeff(:, i))
end
xlabel('frequency (Hz)')
ylabel('PCA coeffs')
legend({'PC1', 'PC2', 'PC3'})
hold off

%%
Features = [mean(PSD(:, 3:16), 2), mean(PSD(:, 17:41), 2), mean(PSD(:, 42:67), 2), mean(PSD(:, 104:300), 2)];
figure
scatter3(Features(:, 1), Features(:, 2), Features(:, 4), 20, colors, 'filled')
%scatter3(score(:, 1), score(:, 2), score(:, 3), 20, colors, 'filled')
xlabel('Delta')
ylabel('Alpha')
zlabel('Gamma')