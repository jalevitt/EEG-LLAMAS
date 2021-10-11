%load('C:\Users\lewislab\Desktop\RealTimeEEGTestData\PhaseTargetingTestData.mat')

highPassFIR = designfilt('highpassfir', 'StopbandFrequency', 0.01 ,'PassbandFrequency', 1, 'SampleRate', EEG.fs);
bandPassFIR = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 50,  'SampleRate', EEG.fs);

WindowsPerSecond = 200;
fs = 200;

TimePoints = {[41:1/WindowsPerSecond:100, 161:1/WindowsPerSecond:220], [101:1/WindowsPerSecond:160, 221:1/WindowsPerSecond:280];
              [111:1/WindowsPerSecond:170, 231:1/WindowsPerSecond:290], [171:1/WindowsPerSecond:230, 291:1/WindowsPerSecond:350]
              [291:1/WindowsPerSecond:320, 351:1/WindowsPerSecond:380], [321:1/WindowsPerSecond:350, 381:1/WindowsPerSecond:410];
}; 
PhaseBins = zeros(1, 19);
PhaseVec = [];
count = 1;
use = [2, 4, 6];
%use = [1, 3, 5];
use = [1];
for i = 1:length(use)
    
    tp = TimePoints{round(i/2), 1};
    data = TestData{use(i)};
    data(:, 1) = filtfilt(bandPassFIR, data(:, 1));
    data(:, 1) = data(:, 1) - mean(data(:, 1));
    phase = angle(hilbert(data(:, 1)));

    idx = 1:length(data(:, 1));
    idx = idx(data(:, 2) == 33022 | data(:, 2) == 33012);

    figure
    plot(data(:, 1));
    hold on
    plot(((data(:, 2) == 33022) | data(:, 2) == 33012) * 20);
    plot(phase * 10);

    for n = 1:length(idx)
        
        if sum(idx(n)/WindowsPerSecond == tp) >= 0 
            p = phase(idx(n));
            PhaseVec(count) = p;
            count = count + 1;
            p = round(9.499 * p/pi) + 10;
            PhaseBins(p) = PhaseBins(p) + 1;
        end
    end

    
    figure
    plot(phase)
end

figure
histogram(PhaseVec, linspace(-1 * pi, pi, 18))
xlabel('Phase (radians)')
ylabel('Hits')
title('Alpha Phase Targeting')
set(gca, 'FontSize', 16)

%%
highPassFIR = designfilt('highpassfir', 'StopbandFrequency', 0.01 ,'PassbandFrequency', 1, 'SampleRate', EEG.fs);
bandPassFIR = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 50,  'SampleRate', EEG.fs);

load('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\Test_12_9_19_base2.mat')
fs = EEG.fs;
idx = EEG.Recording(:, 65) == 33022;
%idx(1:fs*3) = 0;
filtData = filtfilt(bandPassFIR, EEG.Recording(:, EEG.PrimaryChannel));

filtData = filtData - mean(filtData);
phase = angle(hilbert(filtData));

PhaseVec = phase(idx);

figure
histogram(PhaseVec, linspace(-1 * pi, pi, 18))
xlabel('Phase (radians)')
ylabel('Hits')
title('Alpha Phase Targeting')
set(gca, 'FontSize', 16)

TriggerVec = zeros(1000000, 1);
TriggerVec(idx) = 20;
figure
plot(linspace(1/fs, 1000000/fs, 1000000), filtData)
hold on
plot(linspace(1/fs, 1000000/fs, 1000000), TriggerVec)
hold off
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0, 30])

%%
%load('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\Test_11_7_19_sleep1_replay2.mat')

BandPass_SlowWave = designfilt('bandpassiir', ...
    'PassbandFrequency1', 0.5, ...
    'Passbandfrequency2', 2, ...
    'StopbandFrequency1', 0.1 ,...
    'StopbandFrequency2', 10, ...
    'StopbandAttenuation1', 20, ...
    'StopbandAttenuation2', 20, ...
    'PassbandRipple', 0.1, ...
    'DesignMethod', 'butter', ...
    'SampleRate', EEG.fs);
    [b, a] = tf(BandPass_SlowWave);
EEG.Recording = EEG.Recording(1:find(EEG.Recording(:, 1), 1, 'last'), :);
%filt = filtfilt(BandPass_SlowWave, EEG.Recording);
filt = filtfilt(b, a, EEG.Kalman_Signal);
n = size(filt, 1);
ts = (1:n)/EEG.fs;
triggers = find(EEG.Recording(:, end - 1) == 110);
%%
figure
hold on
plot(ts, filt(:, 1))
yline(-1 * EEG.Threshold);
yline(EEG.Threshold);
for i = 1:length(triggers)
    xline(triggers(i)/EEG.fs);
end
xlabel('time (s)')
ylabel('amp (uV)')
xlim([800, 900])
%%
phase = angle(hilbert(filt(:, EEG.PrimaryChannel)));
TriggerPhase = phase(triggers);
figure
histogram(TriggerPhase)
xlabel('Phase (rad)')
ylabel('Count')