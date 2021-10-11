[EEG] = LoadMFF('D:\data\run3_20161007_124807.mff');

%%

EEG.Recording_orig = EEG.Recording;
[n, numChans] = size(EEG.Recording);
EEG.fs_orig = 1000;
EEG.fs = 200;
%%
alphaChannel = 10;

DSrate = EEG.fs_orig/EEG.fs;
dsBuffer = 0;

Copy_orig = zeros(size(EEG.Recording_orig));
EEG.Recording = zeros(ceil(n/5), numChans);


load(strcat('D:\data\run03\hdr.mat'))
bcg_idx = strcmp(hdr.chantype, 'bcg');
eeg_idx = strcmp(hdr.chantype, 'eeg');
temp = 1:length(hdr.chantype);
eeg_chans = temp(eeg_idx);
numEEG = sum(eeg_idx);
KalmanTargets = eeg_chans(randsample(numEEG, 10));
numKalman = length(KalmanTargets);
numBCG = sum(bcg_idx);
numReg = 20;
r = randsample(numBCG, numReg);
x_hat = zeros(numReg + 1, numKalman);
R = 1;
Q = eye(numReg + 1) * 1 * 10 ^ -9;
P = eye(numReg + 1);
P = repmat(P, 1, 1, numKalman);
EEG.Kalman_Signal = zeros(ceil(n/5), numKalman);
EEG.KalmanTargets = KalmanTargets;
EEG.KalmanRegressors = r;

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
DS_FiltSize = round(EEG.fs_orig/20);
DS_Lowpass = firls(DS_FiltSize, [0, 50, 55, EEG.fs_orig/2]/EEG.fs_orig * 2, [1, 1, 0, 0]);
zi = zeros(length(DS_Lowpass) - 1, numChans - 1);

TriggerThreshold = 25;
EEG.Threshold = TriggerThreshold;
LastStimPosition = 0;
z = zeros(6, 1); %filter initial conditions
TriggerBuffer = EEG.fs;
SlowWaveDelay = 0.04;
EEG.StimTimes = zeros(10000, 1);
StimCount = 1;

ntr = 20;
trGap = 0;
        
currentPosition = 1;
currentPosition_orig = 1;
while currentPosition_orig <= n - 100
    ChunkSize = randi([1, 5]) * 10;
    OrigChunk = EEG.Recording_orig(currentPosition_orig:currentPosition_orig + ChunkSize, :)';
    [SamplesInChunk, ChansInChunk] = size(OrigChunk');
    Copy_orig(currentPosition_orig:currentPosition_orig + SamplesInChunk - 1, :) = OrigChunk';
    currentPosition_orig = currentPosition_orig + SamplesInChunk
    
    
    
    
    if trGap~= 0 || sum(Copy_orig(:, ChansInChunk) == 1) > (ntr + 1)
        if trGap == 0
            trSamp = 1:(currentPosition_orig - 1);
            trSamp = trSamp(Copy_orig(:, ChansInChunk) == 1);
            trGap = mode(diff(trSamp))
            EEG.trGap = trGap;
        end
        %template = zeros(trGap, ChansInChunk - 1,  ntr);
        template = reshape(Copy_orig(currentPosition_orig - ntr * trGap:currentPosition_orig - 1, ...
            1:ChansInChunk - 1), trGap, ntr, ChansInChunk - 1);
        template = permute(template, [1, 3, 2]);
        MeanTemplate = mean(template, 3);
        if SamplesInChunk <= trGap
            OrigChunk(1:end - 1, :) = OrigChunk(1:end - 1, :) - ...
                MeanTemplate(end - SamplesInChunk + 1:end, :)';
        else
            numReps = floor(SamplesInChunk / trGap) + 1;
            MeanTemplate = repmat(MeanTemplate, numReps, 1);
            OrigChunk(1:end - 1, :) = OrigChunk(1:end - 1, :) - ...
                MeanTemplate(end - SamplesInChunk + 1:end, :)';

        end
    end
%     [O, zi] = filter(DS_Lowpass, 1, OrigChunk', zi, 1);
%     OrigChunk = O';
      [O, zi] = filter(DS_Lowpass, 1, OrigChunk(1:end - 1, :)', zi, 1);
      OrigChunk(1:end - 1, :) = O';
    
    [chunk, dsBuffer] = DownSampleTriggrs(OrigChunk, DSrate, dsBuffer);
    [SamplesInChunk, ChansInChunk] = size(chunk');
    
    [K_chunk, P, x_hat] = Kalman_chunk(chunk, P, Q, R, r, bcg_idx, x_hat, KalmanTargets);
    EEG.Kalman_Signal(currentPosition:currentPosition + SamplesInChunk - 1, :) ...
        = K_chunk';
    
    EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, :) = chunk';
    %update our position
    currentPosition = currentPosition + SamplesInChunk;
    
    
    sample =  EEG.Kalman_Signal(currentPosition - SamplesInChunk:currentPosition - 1, alphaChannel);
    [FiltSample, z] = filter(b, a, sample', z);

    if max(FiltSample) > TriggerThreshold  && (currentPosition - TriggerBuffer) > LastStimPosition
        %PsychPortAudio('Start', audio_port, repetitions, GetSecs() + SlowWaveDelay, 0);
        %sound(Sound, fsSound)
        EEG.StimTimes(StimCount) = currentPosition;
        StimCount = StimCount + 1;
        LastStimPosition = currentPosition;
    end

   
    
end


%%

n_ds = size(EEG.Recording, 1);
t = (1:n_ds)/EEG.fs;
figure
plot(t, EEG.Recording(:, 10))
xlabel('Time (s)')
ylabel('Voltage (uV)')

figure
plot(t, EEG.Kalman_Signal(:, 10))
xlabel('Time (s)')
ylabel('Voltage (uV)')

%%
EEG.StimTimes = EEG.StimTimes(find(EEG.StimTimes));
filt = filtfilt(b, a, EEG.Kalman_Signal);

phase = angle(hilbert(filt(:, alphaChannel)));
StimPhase = phase(EEG.StimTimes);

figure
histogram(StimPhase)
xlabel('Phase (rads)')
ylabel('Count')

%%
figure
hold on
plot(t, filt(:, alphaChannel))
yline(-1 * EEG.Threshold);
yline(EEG.Threshold);
for i = 1:length(EEG.StimTimes)
    xline(EEG.StimTimes(i)/EEG.fs);
end
xlabel('time (s)')
ylabel('amp (uV)')
xlim([1000, 1100])