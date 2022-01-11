% get raw file
tic
[pathRAW, fileRAW]  = uigetfile({'*.txt'});
[pathEEG, fileEEG]  = uigetfile({'*.mat'});
flicker = input('Enter flicker frequency (Hz): ');
chan = 14;
DS_FiltSize = round(5000 / 20);
DS_Lowpass = firls(DS_FiltSize, [0, 50, 55, 5000 / 2] / 5000 * 2, [1, 1, 0, 0]);
%%
Raw = csvread([fileRAW, pathRAW]);
load([fileEEG, pathEEG]);
toc
RefMatrix = FindRefChans(EEG.KalmanTargets, EEG.KalmanRegressors, 20);
EEG.Recording_orig = Raw;
EEG = OfflineGAC(EEG);
toc
EEG.OfflineGAC = filter(DS_Lowpass, 1, double(EEG.OfflineGAC));
EEG.OfflineDS = downsample(EEG.OfflineGAC, 25);

fs = 200;
StimDur = 16 * fs;
start = find(EEG.Recording(:, end - 1) == EEG.vars.trMarker, 1, 'first');
StimTimes = (0:16) * StimDur + start;

R = 10.^[0:2:10];
Q = 10.^[-12:2:-2];
FiltData = cell(length(R), length(Q));
PowOn = zeros(length(R), length(Q));
PowOff = zeros(length(R), length(Q));
for r = 1:length(R)
    r
    parfor q = 1:length(Q)
        SOff_Pow = 0;
        SOn_Pow = 0;
        K_off = KalmanOffline(EEG.OfflineDS, EEG.KalmanTargets, Q(q), R(r), RefMatrix);
        for i = 2:((length(StimTimes) - 1)/2 - 1)
           pw_off = pwelch(K_off(StimTimes((i - 1) * 2 + 1):StimTimes((i - 1) * 2 + 2), ...
               chan), 200, 100, 1:25, 200);
           pw_on = pwelch(K_off(StimTimes((i - 1) * 2 + 2):StimTimes((i - 1) * 2 + 3), ...
               chan), 200, 100, 1:25, 200);
           SOff_Pow = SOff_Pow + abs(pw_off);
           SOn_Pow = SOn_Pow + abs(pw_on);
        end
        PowOn(r, q) = mag2db(SOn_Pow(flicker) / ((length(StimTimes) - 1)/2));
        PowOff(r, q) = mag2db(SOff_Pow(flicker) / ((length(StimTimes) - 1)/2));
        FiltData{r, q} = K_off;
    end
end
toc
%%
count = 1;
figure
for r = 1:length(R)
    for q = 1:length(Q)
        subplot(length(R), length(Q), count)
        count = count + 1;
        Spec_detrend(FiltData{r, q}(:, chan), 200, 100, 0:0.1:30, fs);
        colormap jet
        caxis([-40, 20])
        title(sprintf('R = %g, Q = %g, Flicker = %d', R(r), Q(q), flicker))
    end
end


%%
count = 1;
figure
for r = 1:length(R)
    for q = 1:length(Q)
        subplot(length(R), length(Q), count)
        count = count + 1;
        bar([1, 2], [PowOn(r, q), PowOff(r, q)])
        xlim([0, 3])
        xticks([1,2])
        xticklabels({'Stim on', 'Stim off'})
        title(sprintf('R = %g, Q = %g, Flicker = %d', R(r), Q(q), flicker))
    end
end
%%
mx = max(max(PowOn - PowOff));
recR = 0;
recQ = 0;
for r = 1:length(R)
    for q = 1:length(Q)
        if PowOn(r, q) - PowOff(r, q) == mx
            recR = r;
            recQ = q;
        end
    end
end
PowOn - PowOff
fprintf('Reccomended R is %g\n', R(recR))
fprintf('Reccomended Q is %g\n', Q(recQ))
toc