[fname, pathname] = uigetfile({'*.mat'});

load([pathname fname]);
display('building filter')
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
display('filtering the signal')
PrimaryChannel = 18;
Signal = EEG.Recording(:, 18);
%%
FiltSignal = Signal;%filter(BandPassFIR, Signal);

time = linspace(1/EEG.fs, length(Signal)/EEG.fs, length(Signal));
display('plotting...')
figure
plot(time, FiltSignal)
xlabel('Time (s)')
ylabel('Filtered Voltage')
pks = findpeaks(FiltSignal);
pks = sort(pks, 'descend');

meanPkHeight = mean(pks(2:10));
display(['The mean peak height is: ', num2str(meanPkHeight)])
display(['The reccomended Threshold is: ', num2str(meanPkHeight * 1.4)])
hold on
line([0, max(time)], [meanPkHeight, meanPkHeight], 'Color', [1, 0, 0])
line([0, max(time)], 1.5 * [meanPkHeight, meanPkHeight], 'Color', [0, 1, 0])
xlim([0, 60 * 5])
ylim([-3 * meanPkHeight, 3 * meanPkHeight])