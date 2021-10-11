% Clear our workspace

if exist('port', 'var')
    fclose(port);
end

clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  DOUBLE CHECK THESE VARIABLES BEFORE YOU START      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG.numChans = 64; %number of channels (not counting triggers)
UseTriggers = 0; %whether or not to use triggers and emit sounds
UseAlphaStim = 0; %whether or not to deploy triggers on alpha waves
UseKalman = 0; %whether or not to use Kalman Filtering
UseSlowWaveStim = 0; %whether or not to deploy triggers on slow waves
fs = 5000; % the native sampling rate
useDownSampling = 1; %determines if we are using downsampling. recommended for new recordings.
targetFs = 200; %the rate to downsample to
alphaChannel = 22; %the primary channel to be used to generating triggers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

%load eyes open/closed classifier
load('EOECDist.mat')
load('EOECModel.mat')

% resolve a stream...
disp('Resolving an EEG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG'); 
end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

% create the channel names
generateLabelNames;

%set up the trigger box
if UseTriggers
    port = serial('COM4');
    fopen(port);
end
IntrinsicLag = 10/200; %This is lag associated with sound production in seconds.

numChans = EEG.numChans;
currentPosition = 1;

%Initiate a blank recording
RecSize = 1000000;
BlockSize = RecSize;
EEG.Recording = zeros(RecSize, numChans + 1);
if UseTriggers
    EEG.Recording = zeros(RecSize, numChans + 2);
end
EEG.fs = fs;
EEG.PrimaryChannel = alphaChannel;
WindowSize = 10;

  % set up Kalman parameters
if UseKalman
    path = uigetdir(); %get path to a hdr file to tell us which channels to use
    load(strcat(path, '/', 'hdr.mat'))
    bcg_idx = strcmp(hdr.chantype, 'bcg');
    eeg_idx = strcmp(hdr.chantype, 'eeg');
    temp = 1:length(hdr.chantype);
    eeg_chans = temp(eeg_idx);
    numEEG = sum(eeg_idx);
    KalmanTargets = eeg_chans(randsample(numEEG, 10));
    numKalman = length(KalmanTargets);
    numBCG = sum(bcg_idx);
    numReg = 30;
    r = randsample(numBCG, numReg);
    x_hat = zeros(numReg + 1, numKalman);
    R = 1;
    Q = eye(numReg + 1) * 1 * 10 ^ -9;
    P = eye(numReg + 1);
    P = repmat(P, 1, 1, numKalman);
    EEG.Kalman_Signal = zeros(RecSize, numKalman);
    EEG.KalmanTargets = KalmanTargets;
    EEG.KalmanRegressors = r;
end

%Set up downsampling
DSrate = EEG.fs/targetFs;
dsBuffer = 0;
if useDownSampling
    EEG.fs = targetFs;
end

%set up filters
highPassIIR = designfilt('highpassiir', 'FilterOrder', 4, 'PassbandFrequency', 5, 'SampleRate', EEG.fs);
highPassFIR = designfilt('highpassfir', 'StopbandFrequency', 0.01 ,'PassbandFrequency', 1, 'SampleRate', EEG.fs);
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

%set up our plots, blank for now
Graph = figure('Position', [10, 50, 1150, 900]);
Eyes = figure();

%set up which channels we will actually plot, by reading from the
%DisplayChannels.csv file
ChannelsToPlot = xlsread('C:\Users\lewislab\Documents\labstreaminglayer\LSL\liblsl-Matlab\RealTimeEEG\DisplayChannels.csv');
if UseTriggers
    ChannelsToPlot = [ChannelsToPlot; 1]; %make sure our trigger channel is on
end
ChannelsToPlot(alphaChannel) = 1; %make sure our primary channel is on
ChannelsToPlot = ChannelsToPlot == 1; %convert to boolean array
numChannelsToPlot = sum(ChannelsToPlot);
EEG.ChannelNames = ChannelNames;

%determine the bin numbers of the frequencies of interest
[psd, bins] = pwelch(rand(EEG.fs,1), [], [], [], EEG.fs);
[~, fullMin] = min(abs(bins - 3));
[~, fullMax] = min(abs(bins - 30));
[~, alphaMin] = min(abs(bins - 8));
[~, alphaMax] = min(abs(bins - 12));

%set up when and how often we'll update our plot
EEGPlotPosition = 0.6 * EEG.fs;
EEGPlotPerSecond = 1;
        
%set up how often we will check for trigger opportunities
alphaPosition = 5 * EEG.fs;
alphaPerSecond = 10;
alphaVec = zeros(1000 * alphaPerSecond, 1);
sawTooth = linspace(0,1,5);
smoothedAlphaVec = zeros(1000 * alphaPerSecond, 1);
pEyesClosed = zeros(1000 * alphaPerSecond, 1);

%build the sound we'll play
fsSound = 48000;
s = 0.05;
frequency = 440;
timeSound = linspace(0, s, fsSound * s + 1);
Sound = hamming(fsSound * s + 1)'.*sin(timeSound * 2 * pi * frequency); %this is just a constant note
cn = dsp.ColoredNoise(1, fsSound * s, 1, 'OutputDataType', 'double');
PinkSound = cn();
PinkSound = hamming(length(PinkSound)).*PinkSound;

%initialize the audio player
InitializePsychSound;
audiopath = 'C:\Users\sdwilli\Desktop\60dbSPL_1000Hz_50mstone_5msramp.wav';
[audiodata, fsSound] = audioread(audiopath);
audio_port = PsychPortAudio('Open', 3, 1, [], fsSound, 2, [], 0.02);
audio_to_play = [PinkSound, PinkSound]';
PsychPortAudio('FillBuffer', audio_port, audio_to_play);
waitForDeviceToStart = 1;
repetitions = 1;

%Initialize storage of our chunk sizes, so we can check latency
EEG.ChunkSizes = zeros(1000000, 1);
count = 1;

%Initialize our slow wave stim variables
TriggerThreshold = 15; 
EEG.Threshold = TriggerThreshold;
LastStimPosition = 0;
z = zeros(6, 1); %filter initial conditions
TriggerBuffer = EEG.fs;
SlowWaveDelay = 0.04;
StimCount = 1;

%Set up Keyboard Reader
[i,j] = GetKeyboardIndices();
KB_task = i(find(strcmp(j, 'Keyboard')));
KbName('UnifyKeyNames');
KbQueueCreate(KB_task);
KbQueueStart(KB_task);
EEG.KeyPresses = zeros(50000, 256);

%now we're ready to start receiving data
disp('Now receiving chunked data...');
EEG.StartTime = GetSecs();
while true
    tic; %start a timer
    % get chunk from the inlet
    [OrigChunk,stamps] = inlet.pull_chunk();
    [SamplesInChunk, ChansInChunk] = size(OrigChunk');

    if numel(OrigChunk) > 0 %if the chunk isn't empty, we need to process it.
        
        %downsample the chunk
        if useDownSampling 
            %down sampling works slightly differently if we're using
            %triggers
            if UseTriggers
                stamps = downsample(stamps', DSrate, dsBuffer)';
                
                %triggers must be downsampled using this specialized
                %function
                [chunk, dsBuffer] = DownSampleTriggrs(OrigChunk, DSrate, dsBuffer);
                
                [SamplesInChunk, ChansInChunk] = size(chunk');
                if SamplesInChunk ~= length(stamps)
                    stamps = zeros(1, SamplesInChunk);
                end
            else
                
                %data without triggers can be downsampled normally
                if SamplesInChunk == 1 && (dsBuffer == 0)
                    chunk = OrigChunk;
                elseif SamplesInChunk == 1
                    chunk = zeros(ChansInChunk, 0);
                else
                    chunk = downsample(OrigChunk', DSrate, dsBuffer)';
                end
                stamps = downsample(stamps', DSrate, dsBuffer)';
                dsBuffer = mod(SamplesInChunk - dsBuffer, DSrate);
                [SamplesInChunk, ~] = size(chunk');
            end
        else
            chunk = OrigChunk;
        end
        %add extra space in our recoding if necesary
        if currentPosition + SamplesInChunk >= RecSize
            EEG.Recording = cat(1, EEG.Recording, zeros(BlockSize, numChans + 1));
            RecSize = RecSize + BlockSize;
            if UseKalman
                 EEG.Kalman_Signal = cat(1, EEG.Kalman_Signal, zeros(BlockSize, numKalman + 1));
            end
        end

        % Kalman filter the chunk
        if UseKalman
            [K_chunk, P, x_hat] = Kalman_chunk(chunk, P, Q, R, r, bcg_idx, x_hat, KalmanTargets);
            EEG.Kalman_Signal(currentPosition:currentPosition + SamplesInChunk - 1, :) ...
                = K_chunk';
        end
        % insert new data into recording
        EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, :) = [chunk', stamps'];
        
        %update our position
        currentPosition = currentPosition + SamplesInChunk;
        
        %keep track of chunk sizes
        EEG.ChunkSizes(count) = SamplesInChunk;
        count = count + 1;
        
        %Read in keyboard presses every 10 chunks
        if mod(count, 10) == 0
            [~, firstpress] = KbQueueCheck(KB_task);
            EEG.KeyPresses(count / 10, :) = firstpress;
        end
        
        %Play a train of sounds at the very beginning just to make sure
        %everything is working
        if sum(count == [10, 30, 50, 70, 90]) > 0
                PsychPortAudio('Start', audio_port, repetitions, GetSecs(), 0);
        end
        
        
        %determine if its tim to update our plot
        if currentPosition > EEGPlotPosition
        %print chunk size once in while so you can monitor it in real time
        %Ideally shoud be in the single digits or low teens but spikes up to ~25 aren't
        %a big deal
            SamplesInChunk
         %calc limits of graph
            xMin = WindowSize * floor((currentPosition - 1)/(WindowSize*EEG.fs));
            xMax = xMin + WindowSize;

            %choose samples to plot
            SampleMin = xMin * EEG.fs + 1;
            SampleMax = currentPosition - 1;

            % make our time vector
            time = linspace(SampleMin/EEG.fs, SampleMax/EEG.fs, SampleMax - SampleMin + 1);
            set(0, 'CurrentFigure', Graph)
            
            %choose the sample section we'll be plotting
            if UseKalman
                SampleToPlot = EEG.Kalman_Signal(SampleMin:SampleMax, :);
            else
                SampleToPlot = EEG.Recording(SampleMin:SampleMax, ChannelsToPlot');
            end
            
            %adjust our sample to create vertical channel offsets
            for i = 1:numChannelsToPlot
                mx = max(SampleToPlot(:, i));
                mn = min(SampleToPlot(:, i));
                SampleToPlot(:, i) = 2 * (SampleToPlot(:, i) - mn)/(mx - mn) + i * 2;
            end
            
            %make our graph

            plot(time, SampleToPlot, 'k' )

            xlim([xMin, xMax])
            ylim([0, numChannelsToPlot * 2 + 2])
            yticks(1 + (2:2:numChannelsToPlot * 2));
            if ~UseKalman
                yticklabels(ChannelNames(ChannelsToPlot))
            end
            xlabel('Time (S)')
            
            %set the next time we'll update our plot
            EEGPlotPosition = EEGPlotPosition + EEG.fs/EEGPlotPerSecond;
        end
        
        if UseSlowWaveStim && SamplesInChunk > 0 
            sample =  EEG.Recording(currentPosition - SamplesInChunk:currentPosition - 1, alphaChannel);
            [FiltSample, z] = filter(b, a, sample', z);

            if max(FiltSample) > TriggerThreshold  && (currentPosition - TriggerBuffer) > LastStimPosition
                PsychPortAudio('Start', audio_port, repetitions, GetSecs() + SlowWaveDelay, 0);
                LastStimPosition = currentPosition;
                EEG.StimTimes(StimCount) = currentPosition;
                StimCount = StimCount + 1;
            end
        end
    end
    
    % check if its time to send new alpha triggers
    if UseAlphaStim && currentPosition > alphaPosition 
        %select our sample and filter it
        sampleForAlpha = EEG.Recording(currentPosition - EEG.fs + 1:currentPosition, alphaChannel);
        FiltData = filtfilt(highPassIIR, sampleForAlpha);
        
        %calc PSD using welch's method
        [psd, ~] = pwelch(FiltData, [], [], [], EEG.fs);
        
        %normalize PSD
        psd = psd/sum(psd(fullMin:fullMax));
        
        %calculate band power, and store it, and a smoothed version
        alphaIDX = alphaPosition/(EEG.fs / alphaPerSecond);
        alphaVec(alphaIDX) = mean(psd(alphaMin:alphaMax));
        smoothedAlphaVec(alphaIDX) = dot(alphaVec(alphaIDX - 4:alphaIDX), sawTooth);
            
        %Determine if the subject's eyes are closed
        [~, score] = predict(EOECModel,  smoothedAlphaVec(alphaIDX));
        pEyesClosed(alphaIDX) = score(1) / sum(score);
        
        %if their eyes are closed, we can look for a chance to emit a sound
        if UseTriggers && pEyesClosed(alphaIDX) > 0.5
            
            %get latest quarter second of data
            FiltData = sampleForAlpha(end - EEG.fs/4:end);
            %filter it
            FiltData = filtfilt(bandPassFIR, FiltData);
            
            %Use the mansouri phase predition algorithm to create a
            %predition 1/2 second long
            [Prediction] = MansouriPhase(FiltData', EEG.fs, EEG.fs / 2);
            
            % read from timer to estimate how much time has passed since
            % chunk was received
            compTime = toc;
            
            % find the delay until the next upcoming peak in our prediction, after ignoring
            % samples lost to compTime and sound prodution lag
            [~, delay] = max(Prediction(round(compTime * EEG.fs + 1 + IntrinsicLag * fs):end));
            
            %convert delay to seconds
            delay = delay/EEG.fs;
            
            %Play sound with delay
            PsychPortAudio('Start', audio_port, repetitions, GetSecs() + delay, 0);

        end
        
        %update our graph of p(eyesClosed)
        set(0, 'CurrentFigure', Eyes)
        plot((1:length(alphaVec))/alphaPerSecond, pEyesClosed)
        xlabel('Time (s)')
        ylabel('Probability eyes closed')
        ylim([-0.5,1.5])
        
        alphaPosition = alphaPosition + EEG.fs / alphaPerSecond;
    end
    pause(0.001)
end

%%
%close stuff
EEG.Recording = EEG.Recording(1:currentPosition, :);
if UseKalman
    EEG.Kalman_Signal = EEG.Kalman_Signal(1:currentPosition, :);
end
EEG.ChunkSizes = EEG.ChunkSizes(1:count);
EEG.StimTimes = EEG.StimTimes(1:StimCount);

disp('transmission completed')
[~, n] = size(EEG.Recording);
EEG.StartTime = EEG.Recording(1, n); %record the start time from the stamps column of our recording
KbQueueRelease(KB_task);
PsychPortAudio('Close')
delete(inlet);
if UseTriggers
    fclose(port);
end
FiltSignal = filter(b, a, EEG.Recording(:, 18));
pks = findpeaks(FiltSignal);
pks = sort(pks, 'descend');

meanPkHeight = mean(pks(2:10));
display(['The mean peak height is: ', num2str(meanPkHeight)])
display(['The reccomended Threshold is: ', num2str(meanPkHeight * 1.4)])