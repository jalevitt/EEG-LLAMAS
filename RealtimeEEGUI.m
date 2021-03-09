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
    PsychPortAudio('Close');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  DOUBLE CHECK THESE VARIABLES BEFORE YOU START      %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    vars.UseTriggers = handles.UseTriggers.Value; %whether or not to use triggers and emit sounds
    vars.UseAlphaStim = handles.UseAlphaStim.Value; %whether or not to deploy triggers on alpha waves
    vars.UseSlowWaveStim = handles.UseSlowWaveStim.Value; %whether or not to deploy triggers on slow waves
    vars.UseGAC = handles.GAC.Value; %whether or not to use Gradient Artifact Detection
    UseOld = handles.UseOld.Value; %wether to replay an old .mat file instead of collecting new data
    fs = str2num(handles.NativeFS.String); % the native sampling rate
    useDownSampling = handles.UseDownsampling.Value; %determines if we are using downsampling. recommended for new recordings.
    targetFs = str2num(handles.TargetFS.String); %the rate to downsample to
    alphaChannel = str2num(handles.PrimaryChannel.String); %the primary channel to be used to generating triggers
    vars.UseKalman = handles.Kalman.Value; % if to use Kalman filter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if were doing GAC, we need downsampling, even if its 1:1
    if vars.UseGAC && ~useDownSampling
        useDownSampling = true;
        targetFs = fs;
    end
    % we also need triggers to do GAC
    if vars.UseGAC && ~vars.UseTriggers
        vars.UseTriggers = true;
    end


    % instantiate the library
    disp('Loading the library...');
    lib = lsl_loadlib();


    %load eyes open/closed classifier
    load('EOECDist.mat')
    load('EOECModel.mat')
    vars.EOECModel = EOECModel;
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
    vars.ChannelNames = ChannelNames;
    
    %set up the trigger box
    if vars.UseTriggers
        port = serial('COM4');
        fopen(port);
    end
    vars.IntrinsicLag = 20/200; %This is lag associated with sound production in seconds.
    
    % open old data if we've configured the recording to replay old data
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
    vars.currentPosition = 1;
    
     
    
%     if vars.UseTriggers
%         disp('Resolving a Markers stream...');
%         result = {};
%         while isempty(result)
%             result = lsl_resolve_byprop(lib,'type','Markers'); 
%         end
%         Mrk_info = result{1}.channel_format();
% 
%         % create a new inlet
%         disp('Opening an inlet...');
%         Mrk_inlet = lsl_inlet(result{1});
% 
%         vals.eventBlockSize = 10000;
%         vals.eventCount = 1;
%         EEG.eventTimes = zeros(vals.eventBlockSize, 1);
%         EEG.eventLbls = cell( vals.eventBlockSize, 1);
%     end
    
    %Initiate a blank recording
    RecSize = 100000;
    BlockSize = RecSize;
    EEG.Recording = zeros(RecSize, numChans + 1);
    if vars.UseTriggers
        EEG.Recording = zeros(RecSize, numChans + 2);
        
    end
    FullWidth = size(EEG.Recording, 2);
    
    EEG.fs = fs;
    EEG.PrimaryChannel = alphaChannel;
    vars.WindowSize = 10;

    %Set up downsampling
    DSrate = EEG.fs/targetFs;
    dsBuffer = 0;
    txtfs = num2str(EEG.fs);
    
    % open .txt files to write raw data and log files
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
    
    %set up downsampled recording
    if useDownSampling
        EEG.fs_orig = EEG.fs;
        EEG.fs = targetFs;
        EEG.Recording_orig = single(zeros(RecSize * DSrate, FullWidth));
        DS_FiltSize = round(EEG.fs_orig/20);
        %[0, 50, 55, EEG.fs_orig/2]/EEG.fs_orig * 2
        DS_Lowpass = firls(DS_FiltSize, [0, 50, 55, EEG.fs_orig/2]/EEG.fs_orig * 2, [1, 1, 0, 0]);
        zi = zeros(length(DS_Lowpass) - 1, numChans);
    else
        EEG.fs_orig = EEG.fs;
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
    vars.highPassIIR = designfilt('highpassiir', 'FilterOrder', 4, 'PassbandFrequency', 5, 'SampleRate', EEG.fs);
    vars.bandPassFIR = designfilt('bandpassfir', 'PassbandFrequency1', 8, 'Passbandfrequency2', 12, 'StopbandFrequency1',6 ,'StopbandFrequency2', 14, 'FilterOrder', 15,  'SampleRate', EEG.fs);
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
    [vars.b, vars.a] = tf(BandPass_SlowWave);
    
    %set up our plots, blank for now
    Graph = figure('Position', [10, 50, 1150, 900]);
    Eyes = figure();

    %set up which channels we will actually plot, by reading from the
    %DisplayChannels.csv file
    vars.ChannelsToPlot = xlsread('C:\Users\lewislab\Documents\labstreaminglayer\LSL\liblsl-Matlab\RealTimeEEG\DisplayChannels.csv');
    if vars.UseTriggers
        vars.ChannelsToPlot = [vars.ChannelsToPlot; 1]; %make sure our trigger channel is on
    end
    vars.ChannelsToPlot(EEG.PrimaryChannel) = 1; %make sure our primary channel is on
    vars.ChannelsToPlot = vars.ChannelsToPlot == 1; %convert to boolean array
    if length(vars.ChannelsToPlot) > EEG.numChans
        vars.ChannelsToPlot = vars.ChannelsToPlot(1:EEG.numChans);
    end
    vars.numChannelsToPlot = sum(vars.ChannelsToPlot);  
    if handles.Kalman.Value
        vars.numChannelsToPlot = numKalman;
    end

    %determine the bin numbers of the frequencies of interest
    [~, bins] = pwelch(rand(EEG.fs,1), [], [], [], EEG.fs);
    [~, vars.fullMin] = min(abs(bins - 3));
    [~, vars.fullMax] = min(abs(bins - 30));
    [~, vars.alphaMin] = min(abs(bins - 8));
    [~, vars.alphaMax] = min(abs(bins - 12));

    %set up when and how often we'll update our plot
    vars.EEGPlotPosition = 0.6 * EEG.fs;
    vars.EEGPlotPerSecond = 1;

    %set up how often we will check for trigger opportunities
    vars.alphaPosition = 5 * EEG.fs;
    vars.alphaPerSecond = 10;
    vars.alphaVec = zeros(1000 * vars.alphaPerSecond, 1);
    vars.sawTooth = linspace(0,1,5);
    vars.smoothedAlphaVec = zeros(1000 * vars.alphaPerSecond, 1);
    vars.pEyesClosed = zeros(1000 * vars.alphaPerSecond, 1);

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
    InitializePsychSound(1);
    vars.audio_port = PsychPortAudio('Open', 5, 1, 3, fsSound, 2, []);
    audio_to_play = [PinkSound'; PinkSound'];
    PsychPortAudio('FillBuffer', vars.audio_port, audio_to_play);
    waitForDeviceToStart = 1;
    vars.repetitions = 1;

    %Initialize storage of our chunk sizes, so we can check latency
    EEG.ChunkSizes = zeros(1000000, 1);
    count = 1;
    
    %Initialize our slow wave stim variables
    TriggerThreshold = str2num(handles.StimThresh.String); 
    EEG.Threshold = TriggerThreshold;
    vars.LastStimPosition = 0;
    vars.z = zeros(6, 1); %filter initial conditions
    vars.TriggerBuffer = EEG.fs;
    vars.SlowWaveDelay = 0.001;
    vars.StimTimes = zeros(10000, 1);
    vars.StimCount = 1;
    
    %Set up Keyboard Reader
    [i,j] = GetKeyboardIndices();
    KB_task = i(find(strcmp(j, 'Keyboard')));
    KbName('UnifyKeyNames');
    KbQueueCreate(KB_task);
    KbQueueStart(KB_task);
    EEG.KeyPresses = zeros(50000, 256);
    
    %Set up file for log and raw data
    FileID = fopen(fname, 'w');
    LogID = fopen(logname, 'w');
    fspec = repmat('%-3.1f,', [1, FullWidth]);
    fspec = [fspec(1:end - 5), '11.7f\r\n'];
    %write to log file
    WriteLogInfo(handles, LogID, txtdate);
    
    %now we're ready to start receiving data
    disp('Now receiving chunked data...');
    vars.currentPosition_orig = 1;
    vars.TrigTime = zeros(5, 1);
    vars.TrigCount = 1;
    AssignTime = 0;
    vars.ntr = 20;
    vars.trGap = 0;
    if UseOld
        Old_position = 1;
    end
    BreakFlag = 0;
    EEG.StartTime = GetSecs();
    vars.startClock = tic;
    
    while (isgraphics(Eyes) && isgraphics(Graph)) && ~BreakFlag
        vars.clock = tic; %start a timer
        % get chunk from the inlet
        t_offset = toc(vars.startClock);
        [vars.OrigChunk,stamps] = inlet.pull_chunk();
        
%         if vars.UseTriggers
%             [mrks,ts] = Mrk_inlet.pull_sample(0.00005);
%         
%             if length(ts) > 0
%                 if length(ts) + vals.eventCount > length(EEG.eventLbls)
%                     EEG.eventTimes = cat(1, EEG.eventTimes, zeros(vals.eventBlockSize, 1));
%                     EEG.eventLbls = cat(1, EEG.eventLbls, cell(vals.eventBlockSize, 1));
%                 end
%                 EEG.eventTimes(vals.eventCount:vals.eventCount + length(ts) - 1) = ts;
%                 EEG.eventLbls(vals.eventCount:vals.eventCount + length(ts) - 1) = mrks;
%                 vals.eventCount = vals.eventCount + length(ts);
%             end
%         end
        [vars.SamplesInChunk, vars.ChansInChunk] = size(vars.OrigChunk');
        %stamps = linspace(t_offset - (length(stamps) - 1) * 1/EEG.fs_orig, t_offset, length(stamps));
        

        if numel(vars.OrigChunk) > 0 %if the chunk isn't empty, we need to process it.
            
            if (UseOld && n_old > Old_position + vars.SamplesInChunk) && (vars.UseGAC || useDownSampling)
                stamps = downsample(stamps, round(5000/EEG.fs_orig), dsBuffer2);
                [vars.OrigChunk, dsBuffer2] = DownSampleTriggrs(vars.OrigChunk, round(5000/EEG.fs_orig), dsBuffer2);
                
                [vars.SamplesInChunk, vars.ChansInChunk] = size(vars.OrigChunk');
                vars.OrigChunk = [EEG_old.Recording(Old_position:Old_position + vars.SamplesInChunk - 1, 1:numChans)'; vars.OrigChunk(end, :)];
                if vars.UseGAC %preserve trs from original
                    vars.OrigChunk(end, find(EEG_old.Recording(Old_position:Old_position + vars.SamplesInChunk - 1, numChans + 1) == 1)) = 1;
                end
                [vars.SamplesInChunk, vars.ChansInChunk] = size(vars.OrigChunk');
                Old_position = Old_position + vars.SamplesInChunk;
            end

            if (UseOld && n_old < Old_position + 5 * vars.SamplesInChunk)
                BreakFlag = 1;
            end
           
            %downsample the chunk
            if useDownSampling 
                    
                
                EEG.Recording_orig(vars.currentPosition_orig:vars.currentPosition_orig + vars.SamplesInChunk - 1, :)...
                    = single([vars.OrigChunk', stamps']);
                
                vars.currentPosition_orig = vars.currentPosition_orig + vars.SamplesInChunk;
                fprintf(FileID, fspec, [vars.OrigChunk', stamps']');
                
                vars = GAC(EEG, vars);

                if vars.UseTriggers
                   [O, zi] = filter(DS_Lowpass, 1, vars.OrigChunk(1:end - 1, :)', zi, 1);
                   vars.OrigChunk(1:end - 1, :) = O';
                else
                    [O, zi] = filter(DS_Lowpass, 1, vars.OrigChunk', zi, 1);
                    vars.OrigChunk = O';
                end
                %down sampling works slightly differently if we're using
                %triggers
                if vars.UseTriggers
                    stamps = downsample(stamps', DSrate, dsBuffer)';

                    %triggers must be downsampled using this specialized
                    %function
                    [chunk, dsBuffer] = DownSampleTriggrs(vars.OrigChunk, DSrate, dsBuffer);

                    [vars.SamplesInChunk, vars.ChansInChunk] = size(chunk');
                    if vars.SamplesInChunk ~= length(stamps)
                        stamps = zeros(1, vars.SamplesInChunk);
                    end
                else

                    %data without triggers can be downsampled normally
                    if vars.SamplesInChunk == 1 && (dsBuffer == 0)
                        chunk = vars.OrigChunk;
                    elseif vars.SamplesInChunk == 1
                        chunk = zeros(vars.ChansInChunk, 0);
                    else
                        chunk = downsample(vars.OrigChunk', DSrate, dsBuffer)';
                    end
                    stamps = downsample(stamps', DSrate, dsBuffer)';
                    dsBuffer = mod(vars.SamplesInChunk - dsBuffer, DSrate);
                    [vars.SamplesInChunk, ~] = size(chunk');
                end
                
            else
                
                chunk = vars.OrigChunk;
                fprintf(FileID, fspec, [chunk', stamps']');
            end
            
            if (UseOld && (n_old > Old_position + vars.SamplesInChunk)) && ~vars.UseGAC
                if EEG.fs < 5000
                    stamps = downsample(stamps', 5000/EEG.fs, dsBuffer)';

                    %triggers must be downsampled using this specialized
                    %function
                    [chunk, dsBuffer] = DownSampleTriggrs(vars.OrigChunk, 5000/EEG.fs, dsBuffer);

                    [vars.SamplesInChunk, vars.ChansInChunk] = size(chunk');
                    if vars.SamplesInChunk ~= length(stamps)
                        stamps = zeros(1, vars.SamplesInChunk);
                    end
                end
                
                chunk = [EEG_old.Recording(Old_position:Old_position + vars.SamplesInChunk - 1, 1:numChans)'; chunk(end, :)];

                Old_position = Old_position + vars.SamplesInChunk;
            end
            %add extra space in our recoding if necesary
            
            if vars.currentPosition + vars.SamplesInChunk >= RecSize
                
                RecSize = RecSize + BlockSize;
                if handles.Kalman.Value
                    EEG.Kalman_Signal = cat(1, EEG.Kalman_Signal, zeros(BlockSize, numKalman));
                end
                if useDownSampling
                    EEG.Recording = cat(1, EEG.Recording, zeros(BlockSize, FullWidth));
                    chop = max(1, vars.currentPosition_orig - EEG.fs_orig * 120);
                    EEG.Recording_orig = cat(1, EEG.Recording_orig(chop:end, :), ...
                        single(zeros(BlockSize * DSrate, FullWidth)));
                    vars.currentPosition_orig = vars.currentPosition_orig - (chop - 1);
                else
                    EEG.Recording = cat(1, EEG.Recording, zeros(BlockSize, FullWidth));
                end
            end
            
            % Kalman filter the chunk
            if handles.Kalman.Value
                [K_chunk, P, x_hat] = Kalman_chunk(chunk, P, Q, R, r, bcg_idx, x_hat, KalmanTargets);
                EEG.Kalman_Signal(vars.currentPosition:vars.currentPosition + vars.SamplesInChunk - 1, :) ...
                    = K_chunk';
            end
            % insert new data into recording
            EEG.Recording(vars.currentPosition:vars.currentPosition + vars.SamplesInChunk - 1, :) = [chunk', stamps'];

            %update our position
            vars.currentPosition = vars.currentPosition + vars.SamplesInChunk;

            %keep track of chunk sizes
            EEG.ChunkSizes(count) = vars.SamplesInChunk;
            count = count + 1;
            
            %periodically send our data to the base workspace so we don't
            %lose it in the event of an error (roughly once a minute)
            %Note that this substantially increases the RAM requirements - 
            %might be usefult to have an option to disable this
            if vars.currentPosition > AssignTime
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
                PsychPortAudio('Start', vars.audio_port, vars.repetitions, GetSecs(), 0);
                vars.TrigTime(vars.TrigCount) = toc(vars.startClock);
                vars.TrigCount = vars.TrigCount  + 1;
            end

            [vars, Graph] = UpdateMainWindow(EEG, vars, Graph);
            [vars] = SlowWaveStim(EEG, vars);
            [vars, Eyes] = AlphaStim(EEG, vars, Eyes);
        end
        
        %This pause seems to be required. I'm not really sure why. You can
        %shorten it if needed, but 1msec is inconsiquential for most
        %purposes
        pause(0.001)
    end
    vars.TrigTime
    disp('Halting recording...')
    fprintf('\n the TRGap is: %d\n', vars.trGap)
    EEG.ChannelNames = vars.ChannelNames;
    
    EEG.Recording = EEG.Recording(1:vars.currentPosition, :);
    if handles.Kalman.Value
        EEG.Kalman_Signal = EEG.Kalman_Signal(1:vars.currentPosition, :);
    end
    EEG.ChunkSizes = EEG.ChunkSizes(1:count);
    EEG.StimTimes = vars.StimTimes(1:vars.StimCount);
    EEG.trGap = vars.trGap;
    
    assignin('base', 'EEG', EEG);
    uisave('EEG')
    
    handles.Recording = EEG;
    guidata(hObject, handles);

    %close stuff
    
    fclose(FileID);
    fclose(LogID);
    if useDownSampling
        %EEG.Recording_orig = EEG.Recording_orig(1:vars.currentPosition_orig, :);
        %EEG.Recording_orig = csvread(fname);
    end
    
    [~, n] = size(EEG.Recording);
    EEG.StartTime = EEG.Recording(1, n); %record the start time from the stamps column of our recording
    KbQueueRelease(KB_task);
    
    PsychPortAudio('Close')
    delete(inlet);
%     delete(Mrk_inlet);
    if vars.UseTriggers
        fclose(port);
    end
    
    
    
    %store our data
    EEG.ChannelNames = vars.ChannelNames;
    assignin('base', 'EEG', EEG);
    assignin('base', 'vars', vars);
    handles.Recording = EEG;
    guidata(hObject, handles);
    disp('transmission completed')
    
    FiltSignal = filter(vars.b, vars.a, EEG.Recording(:, 18));
    if length(FiltSignal) > 1000
        pks = findpeaks(FiltSignal);
        pks = sort(pks, 'descend');

        if length(pks) >= 10
            meanPkHeight = mean(pks(2:10));
            display(['The mean peak height is: ', num2str(meanPkHeight)])
            display(['The reccomended Threshold is: ', num2str(meanPkHeight * 1.4)])
        end
    end



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
