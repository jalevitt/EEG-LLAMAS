
if exist('port', 'var')
    fclose(port);
end


clearvars

load('EOECDist.mat')
load('EOECModel.mat')

% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();


EEG = LoadRealTimeFile('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\Test_11_7_19_sleep1.mat');
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
EEG.Recording = zeros(l, numChans);
Marker_fs = 5000;
WindowSize = 10;

DSrate = Marker_fs/EEG.fs;
dsBuffer = 0;
%%
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
[b, a] = tf(BandPassFIR);
%%
Graph = figure('Position', [10, 100, 1000, 800]);
Eyes = figure();

alphaChannel = 18 + 32;

% data(alphaChannel, :) = zeros(1, l);
% data(alphaChannel, 1:fs:l) = 1;

ChannelsToPlot = rand(numChans - 1, 1) > 0.8;
ChannelsToPlot(numChans - 1) = 1;
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

TriggerThreshold = 15;
LastStimPosition = 0;
z = zeros(6, 1);
TriggerBuffer = EEG.fs;
Delay = 0.04;

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
        EEG.Recording(currentPosition:currentPosition + SamplesInChunk - 1, :) = [data(1:numChans - 2, currentPosition:currentPosition + SamplesInChunk - 1)', chunk(ChansInChunk, :)', stamps'];
        
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
            grid on
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
        if rand < 0.005
            SamplesInChunk
        end
    
        if SamplesInChunk > 0 
            sample =  EEG.Recording(currentPosition - SamplesInChunk:currentPosition - 1, alphaChannel);
            [FiltSample, z] = filter(b, a, sample', z);

            if max(FiltSample) > TriggerThreshold  && (currentPosition - TriggerBuffer) > LastStimPosition
                PsychPortAudio('Start', audio_port, repetitions, GetSecs() + Delay, 0);
                LastStimPosition = currentPosition;
            end
        end
    end
    pause(0.001)
    %PsychPortAudio('Close')
end
%%
disp('transmission completed')
PsychPortAudio('Close')
delete(inlet);
fclose(port);