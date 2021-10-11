function varargout = RealtimeEEGUI(varargin)
% REALTIMEEEGUI MATLAB code for RealtimeEEGUI.fig
%      REALTIMEEEGUI, by itself, creates a new REALTIMEEEGUI or raises the existing
%      singleton*.
%
%      H = REALTIMEEEGUI returns the handle to a new REALTIMEEEGUI or the handle to
%      the existing singleton*.
%
%      REALTIMEEEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REALTIMEEEGUI.M with the given input arguments.
%
%      REALTIMEEEGUI('Property','Value',...) creates a new REALTIMEEEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RealtimeEEGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RealtimeEEGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RealtimeEEGUI

% Last Modified by GUIDE v2.5 17-Sep-2020 10:17:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RealtimeEEGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RealtimeEEGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before RealtimeEEGUI is made visible.
function RealtimeEEGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RealtimeEEGUI (see VARARGIN)

% Choose default command line output for RealtimeEEGUI
handles.output = hObject;
handles.Recording = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RealtimeEEGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RealtimeEEGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Clear our workspace

    if exist('port', 'var')
        fclose(port);
    end
    
    if ~isempty(instrfind)
        fclose(instrfind);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  DOUBLE CHECK THESE VARIABLES BEFORE YOU START      %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    UseTriggers = handles.UseTriggers.Value; %whether or not to use triggers and emit sounds
    UseAlphaStim = handles.UseAlphaStim.Value; %whether or not to deploy triggers on alpha waves
    UseSlowWaveStim = handles.UseSlowWaveStim.Value; %whether or not to deploy triggers on slow waves
    UseGAC = handles.GAC.Value; %whether or not to use Gradient Artifact Detection
    UseOld = handles.UseOld.Value; %wether to replay an old .mat file instead of collecting new data
    fs = str2num(handles.NativeFS.String); % the native sampling rate
    useDownSampling = handles.UseDownsampling.Value; %determines if we are using downsampling. recommended for new recordings.
    targetFs = str2num(handles.TargetFS.String); %the rate to downsample to
    alphaChannel = str2num(handles.PrimaryChannel.String); %the primary channel to be used to generating triggers

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
    
    if UseOld
        
        LoadMethod = menu('How would you like to load a file?' , 'From File', 'From .mff directory');
        if LoadMethod == 1
            [filename, pathname] = uigetfile({'*.mat'});
            load([pathname, filename]);
            handles.ReplayPath = [pathname, filename];
            EEG_old = EEG;
        else
            pathname = uigetdir();
            handles.ReplayPath = pathname;
            EEG_old = LoadMFF(pathname);
        end
        n_old = size(EEG_old.Recording, 1);
        EEG = struct();
        EEG.numChans = size(EEG_old.Recording, 2) - 1;
        dsBuffer2 = 0;
    else
        EEG = struct();
        EEG.numChans = str2num(handles.numChans.String); %number of channels (not counting triggers)
    end
    
    numChans = EEG.numChans;
    currentPosition = 1;
    

    %Initiate a blank recording
    RecSize = 100000;
    BlockSize = RecSize;
    EEG.Recording = zeros(RecSize, numChans + 1);
    if UseTriggers
        EEG.Recording = zeros(RecSize, numChans + 2);
        
    end
    FullWidth = size(EEG.Recording, 2);
    
    EEG.fs = fs;
    EEG.PrimaryChannel = alphaChannel;
    WindowSize = 10;

    %Set up downsampling
    DSrate = EEG.fs/targetFs;
    dsBuffer = 0;
    txtfs = num2str(EEG.fs);
    
    txtdate = datestr(datetime('now'));
    txtdate = replace(txtdate, ':', '_');
    txtdate = replace(txtdate, ' ', '_');
    txtdate = replace(txtdate, '-', '_');
    if isfield(handles, 'StreamPath')
        fname = [handles.StreamPath, txtdate, '_fs_', txtfs, '.txt'];
        logname = [handles.StreamPath, 'LOG', txtdate, '_fs_', txtfs, '.txt'];
    else
        fname = ['C:\Users\lewislab\Documents\labstreaminglayer\LSL\liblsl-Matlab\RealTimeEEG\Recordings\', txtdate, '_fs_', txtfs, '.txt'];
        logname = ['C:\Users\lewislab\Documents\labstreaminglayer\LSL\liblsl-Matlab\RealTimeEEG\Recordings\LOG', txtdate, '_fs_', txtfs, '.txt'];

    end
    
    if useDownSampling
        EEG.fs_orig = EEG.fs;
        EEG.fs = targetFs;
        EEG.Recording_orig = single(zeros(RecSize * DSrate, FullWidth));
        DS_FiltSize = round(EEG.fs_orig/20);
        %[0, 50, 55, EEG.fs_orig/2]/EEG.fs_orig * 2
        DS_Lowpass = firls(DS_FiltSize, [0, 50, 55, EEG.fs_orig/2]/EEG.fs_orig * 2, [1, 1, 0, 0]);
        zi = zeros(length(DS_Lowpass) - 1, numChans);
    end
    
    % set up Kalman parameters
    if handles.Kalman.Value
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
    
    
    %set up filters
    highPassIIR = designfilt('highpassiir', 'FilterOrder', 4, 'PassbandFrequency', 5, 'SampleRate', EEG.fs);
    highPassFIR = designfilt('highpassfir', 'StopbandFrequency', 0.01 ,'PassbandFrequency', 1, 'SampleRate', EEG.fs);
    bandPassFIR = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 15,  'SampleRate', EEG.fs);
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
    if handles.Kalman.Value
        numChannelsToPlot = numKalman;
    end

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
    Sound = hamming(fsSound * s + 1)'.*sin(timeSound * 2 * pi * frequency);
    cn = dsp.ColoredNoise(1, fsSound * s, 1, 'OutputDataType', 'double');
    PinkSound = cn();
    PinkSound = hamming(length(PinkSound)).*PinkSound;

    %initialize the audio player
    InitializePsychSound;
    audio_port = PsychPortAudio('Open', 3, 1, [], fsSound, 2, [], 0.02);
    audio_to_play = [PinkSound'; PinkSound'];
    PsychPortAudio('FillBuffer', audio_port, audio_to_play);
    waitForDeviceToStart = 1;
    repetitions = 1;

    %Initialize storage of our chunk sizes, so we can check latency
    EEG.ChunkSizes = zeros(1000000, 1);
    count = 1;
    
    %Initialize our slow wave stim variables
    TriggerThreshold = str2num(handles.StimThresh.String); 
    EEG.Threshold = TriggerThreshold;
    LastStimPosition = 0;
    z = zeros(6, 1); %filter initial conditions
    TriggerBuffer = EEG.fs;
    SlowWaveDelay = 0.04;
    EEG.StimTimes = zeros(10000, 1);
    StimCount = 1;
    
    %Set up Keyboard Reader
    [i,j] = GetKeyboardIndices();
    KB_task = i(find(strcmp(j, 'Keyboard')));
    KbName('UnifyKeyNames');
    KbQueueCreate(KB_task);
    KbQueueStart(KB_task);
    EEG.KeyPresses = zeros(50000, 256);
    
    
    FileID = fopen(fname, 'w');
    LogID = fopen(logname, 'w');
    fspec = repmat('%-3.1f,', [1, FullWidth]);
    fspec = [fspec(1:end - 5), '11.7f\r\n'];
    
    WriteLogInfo(handles, LogID, txtdate);
    
    %now we're ready to start receiving data
    disp('Now receiving chunked data...');
    EEG.StartTime = GetSecs();
    currentPosition_orig = 1;
    AssignTime = 0;
    ntr = 20;
    trGap = 0;
    if UseOld
        Old_position = 1;
    end
    BreakFlag = 0;
    while (isgraphics(Eyes) && isgraphics(Graph)) && ~BreakFlag
        tic; %start a timer
        % get chunk from the inlet
        [OrigChunk,stamps] = inlet.pull_chunk();
        [SamplesInChunk, ChansInChunk] = size(OrigChunk');
        

        if numel(OrigChunk) > 0 %if the chunk isn't empty, we need to process it.
            
            if (UseOld && n_old > Old_position + SamplesInChunk) && (UseGAC || useDownSampling)
                stamps = downsample(stamps, round(5000/EEG.fs_orig), dsBuffer2);
                [OrigChunk, dsBuffer2] = DownSampleTriggrs(OrigChunk, round(5000/EEG.fs_orig), dsBuffer2);
                
                [SamplesInChunk, ChansInChunk] = size(OrigChunk');
                OrigChunk = [EEG_old.Recording(Old_position:Old_position + SamplesInChunk - 1, 1:numChans)'; OrigChunk(end, :)];
                if UseGAC %preserve trs from original
                    OrigChunk(end, find(EEG_old.Recording(Old_position:Old_position + SamplesInChunk - 1, numChans + 1) == 1)) = 1;
                end
                [SamplesInChunk, ChansInChunk] = size(OrigChunk');
                Old_position = Old_position + SamplesInChunk;
            end

            if (UseOld && n_old < Old_position + 5 * SamplesInChunk)
                BreakFlag = 1;
            end
           
            %downsample the chunk
            if useDownSampling 
                    
                
                EEG.Recording_orig(currentPosition_orig:currentPosition_orig + SamplesInChunk - 1, :)...
                    = single([OrigChunk', stamps']);
                
                currentPosition_orig = currentPosition_orig + SamplesInChunk;
                fprintf(FileID, fspec, [OrigChunk', stamps']');
                
                if UseGAC % Perform GAC
                    if trGap~= 0 || sum(EEG.Recording_orig(:, end - 1) == 1) > (ntr + 1)
                        if trGap == 0
                            trSamp = 1:(currentPosition_orig - 1);
                            trSamp = trSamp(EEG.Recording_orig(:, end - 1) == 1);
                            trGap = mode(diff(trSamp))
                            EEG.trGap = trGap;
                        end
                        %template = zeros(trGap, ChansInChunk - 1,  ntr);
                        template = reshape(EEG.Recording_orig(currentPosition_orig - ntr * trGap:currentPosition_orig - 1, ...
                            1:end - 2), trGap, ntr, ChansInChunk - 1);
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
                end

                if UseTriggers
                   [O, zi] = filter(DS_Lowpass, 1, OrigChunk(1:end - 1, :)', zi, 1);
                   OrigChunk(1:end - 1, :) = O';
                else
                    [O, zi] = filter(DS_Lowpass, 1, OrigChunk', zi, 1);
                    OrigChunk = O';
                end
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
                fprintf(FileID, fspec, [chunk', stamps']');
            end
            
            if (UseOld && n_old > Old_position + SamplesInChunk) && ~UseGAC
                chunk = [EEG_old.Recording(Old_position:Old_position + SamplesInChunk - 1, 1:numChans)'; chunk(end, :)];

                Old_position = Old_position + SamplesInChunk;
            end
            %add extra space in our recoding if necesary
            
            if currentPosition + SamplesInChunk >= RecSize
                
                RecSize = RecSize + BlockSize;
                if handles.Kalman.Value
                    EEG.Kalman_Signal = cat(1, EEG.Kalman_Signal, zeros(BlockSize, numKalman));
                end
                if useDownSampling
                    EEG.Recording = cat(1, EEG.Recording, zeros(BlockSize, FullWidth));
                    chop = currentPosition_orig - EEG.fs_orig * 120;
                    EEG.Recording_orig = cat(1, EEG.Recording_orig(chop:end, :), ...
                        single(zeros(BlockSize * DSrate, FullWidth)));
                    currentPosition_orig = currentPosition_orig - (chop - 1);
                else
                    EEG.Recording = cat(1, EEG.Recording, zeros(BlockSize, FullWidth));
                end
            end
            
            % Kalman filter the chunk
            if handles.Kalman.Value
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
            
            %periodically send our data to the base workspace so we don't
            %lose it in the event of an error (roughly once a minute)
            if currentPosition > AssignTime
                assignin('base', 'EEG', EEG);
                AssignTime = AssignTime + 60 * EEG.fs;
            end
            
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

            %determine if its time to update our plot
            if currentPosition > EEGPlotPosition
                %print chunk size once in while so you can monitor it in real time
                %Ideally shoud be in the single digits or low teens but
                %spikes up to ~25 aren't unusual
                
                SamplesInChunk
                %calc limits of graph
                
                xMin = WindowSize * floor((currentPosition - 1)/(WindowSize*EEG.fs));
                xMax = xMin + WindowSize;

                %choose samples to plot
                SampleMin = xMin * EEG.fs + 1;
                SampleMax = currentPosition - 1;

                % make our time vector
                time = linspace(SampleMin/EEG.fs, SampleMax/EEG.fs, SampleMax - SampleMin + 1);
                

                %choose the sample section we'll be plotting
                if handles.Kalman.Value
                    SampleToPlot = EEG.Kalman_Signal(SampleMin:SampleMax, :);
                    if UseTriggers
                        SampleToPlot = [SampleToPlot, EEG.Recording(SampleMin:SampleMax, end - 1)];
                    end
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
                if isgraphics(Graph)
                    set(0, 'CurrentFigure', Graph)
                    plot(time, SampleToPlot, 'k' )
                    xLines = unique(round(xMin:xMax));
                    hold on
                    for xl = xLines
                        xline(xl, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.25);
                    end
                    hold off
                    xlim([xMin, xMax])
                    ylim([1, numChannelsToPlot * 2 + 2])
                    yticks(1 + (2:2:numChannelsToPlot * 2));
                    if ~handles.Kalman.Value
                        yticklabels(ChannelNames(ChannelsToPlot))
                    end
                    xlabel('Time (S)')
                end

                %set the next time we'll update our plot
                EEGPlotPosition = EEGPlotPosition + EEG.fs/EEGPlotPerSecond;
            end
            if UseSlowWaveStim && SamplesInChunk > 0 
                if ~handles.Kalman.Value
                    sample =  EEG.Recording(currentPosition - SamplesInChunk:currentPosition - 1, alphaChannel);
                else
                    sample =  EEG.Kalman_Signal(currentPosition - SamplesInChunk:currentPosition - 1, alphaChannel);
                end
                [FiltSample, z] = filter(b, a, sample', z);

                if max(FiltSample) > TriggerThreshold  && (currentPosition - TriggerBuffer) > LastStimPosition
                    PsychPortAudio('Start', audio_port, repetitions, GetSecs() + SlowWaveDelay, 0);
                    %sound(Sound, fsSound)
                    EEG.StimTimes(StimCount) = currentPosition;
                    StimCount = StimCount + 1;
                    LastStimPosition = currentPosition;
                end
            end
        end

        % check if its time to send new triggers
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
                [~, delay] = max(Prediction(round(compTime * EEG.fs + 1 + IntrinsicLag * EEG.fs):end));

                %convert delay to seconds
                delay = delay/EEG.fs;

                %Play sound with delay
                PsychPortAudio('Start', audio_port, repetitions, GetSecs() + delay, 0);
                t = round(currentPosition + IntrinsicLag * EEG.fs + compTime * EEG.fs + delay * EEG.fs);
                if ~numel(t) == 0
                    
                    EEG.StimTimes(StimCount) = t(1);
                else
                    EEG.StimTimes(StimCount) = currentPosition;
                end
                StimCount = StimCount + 1;
                LastStimPosition = currentPosition;

            end

            %update our graph of p(eyesClosed)
            if isgraphics(Eyes)
                set(0, 'CurrentFigure', Eyes)
                plot((1:length(alphaVec))/alphaPerSecond, pEyesClosed)
                hold on
                xline(currentPosition / EEG.fs, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.25);
                hold off
                xlabel('Time (s)')
                ylabel('Probability eyes closed')
                ylim([-0.5,1.5])
            end

            alphaPosition = alphaPosition + EEG.fs / alphaPerSecond;
        end
        pause(0.001)
    end
    
    disp('Halting recording...')
    fprintf('\n the TRGap is: %d\n', trGap)
    EEG.ChannelNames = ChannelNames;
    
    EEG.Recording = EEG.Recording(1:currentPosition, :);
    if handles.Kalman.Value
        EEG.Kalman_Signal = EEG.Kalman_Signal(1:currentPosition, :);
    end
    EEG.ChunkSizes = EEG.ChunkSizes(1:count);
    EEG.StimTimes = EEG.StimTimes(1:StimCount);
    
    assignin('base', 'EEG', EEG);
    uisave('EEG')
    
    handles.Recording = EEG;
    guidata(hObject, handles);

    %close stuff
    
    fclose(FileID);
    fclose(LogID);
    if useDownSampling
        %EEG.Recording_orig = EEG.Recording_orig(1:currentPosition_orig, :);
        %EEG.Recording_orig = csvread(fname);
    end
    
    [~, n] = size(EEG.Recording);
    EEG.StartTime = EEG.Recording(1, n); %record the start time from the stamps column of our recording
    KbQueueRelease(KB_task);
    
    PsychPortAudio('Close')
    delete(inlet);
    if UseTriggers
        fclose(port);
    end
    
    
    
    %store our data
    EEG.ChannelNames = ChannelNames;
    assignin('base', 'EEG', EEG);
    handles.Recording = EEG;
    guidata(hObject, handles);
    disp('transmission completed')
    
    FiltSignal = filter(b, a, EEG.Recording(:, 18));
    pks = findpeaks(FiltSignal);
    pks = sort(pks, 'descend');

    meanPkHeight = mean(pks(2:10));
    display(['The mean peak height is: ', num2str(meanPkHeight)])
    display(['The reccomended Threshold is: ', num2str(meanPkHeight * 1.4)])



% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    EEG = handles.Recording;
    uisave('EEG');
    display('Save Completed!')



function numChans_Callback(hObject, eventdata, handles)
% hObject    handle to numChans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numChans as text
%        str2double(get(hObject,'String')) returns contents of numChans as a double


% --- Executes during object creation, after setting all properties.
function numChans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numChans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NativeFS_Callback(hObject, eventdata, handles)
% hObject    handle to NativeFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NativeFS as text
%        str2double(get(hObject,'String')) returns contents of NativeFS as a double


% --- Executes during object creation, after setting all properties.
function NativeFS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NativeFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TargetFS_Callback(hObject, eventdata, handles)
% hObject    handle to TargetFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetFS as text
%        str2double(get(hObject,'String')) returns contents of TargetFS as a double


% --- Executes during object creation, after setting all properties.
function TargetFS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UseDownsampling.
function UseDownsampling_Callback(hObject, eventdata, handles)
% hObject    handle to UseDownsampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseDownsampling



function PrimaryChannel_Callback(hObject, eventdata, handles)
% hObject    handle to PrimaryChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PrimaryChannel as text
%        str2double(get(hObject,'String')) returns contents of PrimaryChannel as a double


% --- Executes during object creation, after setting all properties.
function PrimaryChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PrimaryChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UseTriggers.
function UseTriggers_Callback(hObject, eventdata, handles)
% hObject    handle to UseTriggers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseTriggers


% --- Executes on button press in UseSlowWaveStim.
function UseSlowWaveStim_Callback(hObject, eventdata, handles)
% hObject    handle to UseSlowWaveStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseSlowWaveStim


% --- Executes on button press in UseAlphaStim.
function UseAlphaStim_Callback(hObject, eventdata, handles)
% hObject    handle to UseAlphaStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseAlphaStim


% --- Executes on button press in Kalman.
function Kalman_Callback(hObject, eventdata, handles)
% hObject    handle to Kalman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Kalman


% --- Executes on button press in AveRef.
function AveRef_Callback(hObject, eventdata, handles)
% hObject    handle to AveRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AveRef


% --- Executes on button press in GAC.
function GAC_Callback(hObject, eventdata, handles)
% hObject    handle to GAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GAC


% --- Executes on button press in UseOld.
function UseOld_Callback(hObject, eventdata, handles)
% hObject    handle to UseOld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseOld


% --- Executes on button press in SetPath.
function SetPath_Callback(hObject, eventdata, handles)
% hObject    handle to SetPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.StreamPath = uigetdir();
fprintf('Set output stream to %s\n', handles.StreamPath)
guidata(hObject, handles);




function StimThresh_Callback(hObject, eventdata, handles)
% hObject    handle to StimThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimThresh as text
%        str2double(get(hObject,'String')) returns contents of StimThresh as a double


% --- Executes during object creation, after setting all properties.
function StimThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
