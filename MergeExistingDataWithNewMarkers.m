
if exist('port', 'var')
    fclose(port);
end


clearvars

load('EOECDist.mat')
load('EOECModel.mat')

% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();


EEG = LoadRealTimeFile('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19.mat');
EEG.data = EEG.data(:, 1: find(EEG.data(1, :), 1, 'last'));
chunkSize = 32;
fs = EEG.srate;
EEG.fs = fs;

data = double(EEG.data);
EEG.data = [];
EEG.Recording = [];
[numChans, l] = size(data);

IntrinsicLag = 10/200;

% resolve a stream...
disp('Resolving an EEG stream...');
result = {};

while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG'); 
end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

generateLabelNames;

EEG.numChans = numChans;

UseTriggers = 1;
if UseTriggers
    port = serial('COM4');
    fopen(port);
end

currentPosition = 1;

RecSize = 1000000;
BlockSize = RecSize;
EEG.Recording = zeros(l, numChans + 1);
Marker_fs = 5000;
WindowSize = 10;

DSrate = Marker_fs/EEG.fs;
dsBuffer = 0;

highPassIIR = designfilt('highpassiir', 'FilterOrder', 4, 'PassbandFrequency', 5, 'SampleRate', EEG.fs);
bandPassFIR = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 15,  'SampleRate', EEG.fs);

Graph = figure('Position', [10, 100, 1000, 800]);
Eyes = figure();

alphaChannel = 22;

% data(alphaChannel, :) = zeros(1, l);
% data(alphaChannel, 1:fs:l) = 1;

ChannelsToPlot = rand(numChans, 1) > 0.8;
ChannelsToPlot(numChans) = 1;
ChannelsToPlot(alphaChannel) = 1;

numChannelsToPlot = sum(ChannelsToPlot);

[psd, bins] = pwelch(rand(EEG.fs,1), [], [], [], EEG.fs);
[~, fullMin] = min(abs(bins - 3));
[~, fullMax] = min(abs(bins - 30));
[~, alphaMin] = min(abs(bins - 8));
[~, alphaMax] = min(abs(bins - 12));

EEGPlotPosition = 0.6 * EEG.fs;
EEGPlotPerSecond = 1;

alphaPosition = 5 * EEG.fs;
alphaPerSecond = 10;
alphaVec = zeros(1000 * alphaPerSecond, 1);
sawTooth = linspace(0,1,5);
smoothedAlphaVec = zeros(1000 * alphaPerSecond, 1);
pEyesClosed = zeros(1000 * alphaPerSecond, 1);

tp = zeros(1000, 1);
count = 1;

fsSound = 48000;
s = 0.05;
frequency = 440;
timeSound = linspace(0, s, fsSound * s + 1);
Sound = hamming(fsSound * s + 1)'.*sin(timeSound * 2 * pi * frequency);
%Sound = sin(timeSound * 2 * pi * frequency);

InitializePsychSound;
audio_port = PsychPortAudio('Open', [], 1, [], fsSound, 2, [], 0.015);
audio_to_play = [Sound; Sound];
PsychPortAudio('FillBuffer', audio_port, audio_to_play);
waitForDeviceToStart = 1;
repetitions = 1;
PsychPortAudio('Start', audio_port, repetitions, 0, waitForDeviceToStart);

disp('Now receiving chunked data...');
while currentPosition < l - fs
        % get chunk from the inlet
    tic;
    [chunk,stamps] = inlet.pull_chunk();
    [SamplesInChunk, ChansInChunk] = size(chunk');
    
    if numel(chunk) > 0
        % downsample the marker chunk
        stamps = downsample(stamps', DSrate, dsBuffer)';
        [chunk, dsBuffer] = DownSampleTriggrs(chunk, DSrate, dsBuffer);
        [SamplesInChunk, ChansInChunk] = size(chunk');
        
        if SamplesInChunk ~= length(stamps)
            stamps = zeros(1, SamplesInChunk);
        end
        
        %plug our chunk into our recording
%         EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, numChans) = chunk(ChansInChunk, :)';
%         EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, 1:numChans - 1) = data(1:numChans - 1, currentPosition:currentPosition + SamplesInChunk - 1)';
%         EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, numChans + 1) = stamps';
        EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, :) = [data(1:numChans - 1, currentPosition:currentPosition + SamplesInChunk - 1)', chunk(ChansInChunk, :)', stamps'];
        
        currentPosition = currentPosition + SamplesInChunk;
        
        if currentPosition > EEGPlotPosition
        %calc limits of graph
            xMin = WindowSize * floor((currentPosition - 1)/(WindowSize*EEG.fs));
            xMax = xMin + WindowSize;

            %choose saples to plot
            SampleMin = xMin * EEG.fs + 1;
            SampleMax = currentPosition - 1;

            %make our graph
            time = linspace(SampleMin/EEG.fs, SampleMax/EEG.fs, SampleMax - SampleMin + 1);
            set(0, 'CurrentFigure', Graph)

            SampleToPlot = EEG.Recording(SampleMin:SampleMax, ChannelsToPlot);
            for i = 1:numChannelsToPlot

                mx = max(SampleToPlot(:, i));
                mn = min(SampleToPlot(:, i));
                SampleToPlot(:, i) = 2 * (SampleToPlot(:, i) - mn)/(mx - mn) + i * 2;

            end
            plot(time, SampleToPlot, 'k' )

            xlim([xMin, xMax])
            ylim([0, numChannelsToPlot * 2 + 2])
            yticks(2:2:numChannelsToPlot * 2);
            yticklabels(ChannelNames(ChannelsToPlot))
            xlabel('Time (S)')
            
            EEGPlotPosition = EEGPlotPosition + EEG.fs/EEGPlotPerSecond;
        end
%          if rand < 0.01
%              soundTemp = audio_to_play;
%              
%              %sound(Sound, fsSound, 8);
%              PsychPortAudio('FillBuffer', audio_port, soundTemp);
%              PsychPortAudio('Start', audio_port, repetitions, 0, waitForDeviceToStart);
%              %fprintf(port, 255);
%              tp(count) = currentPosition + toc * fs;
%              count = count + 1;
%          end
            
    
    end
    if currentPosition > alphaPosition
        sampleForAlpha = EEG.Recording(currentPosition - EEG.fs + 1:currentPosition, alphaChannel);
        FiltData = filtfilt(highPassIIR, sampleForAlpha);
        %FiltData = (FiltData - mean(FiltData)) / std(FiltData);
        %[psd, bins] = periodogram(FiltData, [], [], fs);
        [psd, ~] = pwelch(FiltData, [], [], [], EEG.fs);

        psd = psd/sum(psd(fullMin:fullMax));
        alphaIDX = alphaPosition/(EEG.fs / alphaPerSecond);
        alphaVec(alphaIDX) = mean(psd(alphaMin:alphaMax));
        smoothedAlphaVec(alphaIDX) = dot(alphaVec(alphaIDX - 4:alphaIDX), sawTooth);
        
        [~, score] = predict(EOECModel,  smoothedAlphaVec(alphaIDX));
        pEyesClosed(alphaIDX) = score(1) / sum(score);
        if pEyesClosed(alphaIDX) > 0.5 
%             [~, target] = max(sampleForAlpha);
%             compTime = toc;
%             delay = (target) / fs - compTime - IntrinsicLag;
% 
%             if delay <= 0
%                 delay = delay + 1;
%             end
            FiltData = sampleForAlpha(end - EEG.fs/4:end);
            FiltData = filtfilt(bandPassFIR, FiltData);
            [Prediction] = MansouriPhase(FiltData', EEG.fs, EEG.fs / 2);
            compTime = toc;
            [~, delay] = min(Prediction(round(compTime * EEG.fs + 1 + IntrinsicLag * fs):end));
            
            delay = delay/EEG.fs;
            
            %soundTemp = [zeros(2, round(delay * fsSound)), audio_to_play];
            
%             audio_port = PsychPortAudio('Open', [], 1, [], fsSound, 2, [], delay);
%             PsychPortAudio('FillBuffer', audio_port, audio_to_play);
            PsychPortAudio('Start', audio_port, repetitions, GetSecs() + delay, 0);%waitForDeviceToStart);
            %sound([zeros(1, round(delay * fsSound)), Sound], fsSound);
            %sound(Sound, fsSound);
            %fprintf(port, 255);
        end

        
        set(0, 'CurrentFigure', Eyes)
        plot((1:length(alphaVec))/alphaPerSecond, pEyesClosed)
        xlabel('Time (s)')
        ylabel('Probability eyes closed')
        ylim([-0.5,1.5])
        
        alphaPosition = alphaPosition + EEG.fs / alphaPerSecond;
    end
    pause(0.001)
    %PsychPortAudio('Close')
end

disp('transmission completed')
PsychPortAudio('Close')
delete(inlet);
fclose(port);