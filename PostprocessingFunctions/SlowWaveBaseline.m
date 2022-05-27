function [vars, Graph, EEG] = SlowWaveBaseline(EEG, vars, Graph)
% obtain baseline measures for slow wave stim

if vars.SamplesInChunk > 0
    if ~isfield(vars, 'PhasePredictor')
        if EEG.PrimaryChannel == 17
            fprintf('loading Fpz predictor...')
            load('03-08-2022 15-35results_Fpz_2subs_100delay.mat', 'results');
        elseif EEG.PrimaryChannel == 6
            fprintf('loading C4 predictor...')
            load('04-20-2022 12-33_C4_3subs_100ms.mat', 'results');
        else
            fprintf('WARNING! NO PREDICTOR FOR THIS CHANNEL')
        end
        vars.PhasePredictor = resetState(results(1).net);
        vars.SlowWaveDelay = .000;
        vars.Angles = zeros(1000000, 1);
        vars.X = zeros(3, 1000000);
    end
    if ~isfield(vars, 'b_delta')
        
        vars.delpthresh = 7; % must be higher than this value
        vars.movsthresh = 4; % must be lower than this value
        vars.TriggerBuffer = EEG.fs.*2;
        
        vars.allMags = 0;
        vars.alldelps = 0;
        vars.allmovs = 0;

        BandPass_SlowWave = designfilt('bandpassiir', ...
            'PassbandFrequency1', 0.4, ... 
            'Passbandfrequency2', 1.2, ... 
            'StopbandFrequency1', 0.01,...
            'StopbandFrequency2', 4, ...
            'StopbandAttenuation1', 30, ...
            'StopbandAttenuation2', 30, ...
            'PassbandRipple', 1, ...
            'DesignMethod', 'butter', ...
            'SampleRate', EEG.fs);
        [vars.b, vars.a] = tf(BandPass_SlowWave);
        
        delta_filter = designfilt('bandpassiir', ...
            'PassbandFrequency1', 1, ... 
            'Passbandfrequency2', 4, ... 
            'StopbandFrequency1', .01,...
            'StopbandFrequency2', 9, ...
            'StopbandAttenuation1', 20, ...
            'StopbandAttenuation2', 20, ...
            'PassbandRipple', 1, ...
            'DesignMethod', 'butter', ...
            'SampleRate', EEG.fs);
        [vars.b_delta, vars.a_delta] = tf(delta_filter);
        
        mov_filter = designfilt('highpassiir', ...
            'PassbandFrequency',20, ... 
            'StopbandFrequency', 5,...
            'StopbandAttenuation', 30, ...
            'PassbandRipple', 1, ...
            'DesignMethod', 'butter', ...
            'SampleRate', EEG.fs);
        [vars.b_mov, vars.a_mov] = tf(mov_filter);
        
        hp_filter = designfilt('highpassiir', ...
            'PassbandFrequency', 0.4, ... 
            'StopbandFrequency', 0.01,...
            'StopbandAttenuation', 40, ...
            'PassbandRipple', 0.1, ...
            'DesignMethod', 'butter', ...
            'SampleRate', EEG.fs);
        [vars.b_hp, vars.a_hp] = tf(hp_filter);
        vars.zhp = zeros(2,1); %filter initial conditions  
    end
    if ~vars.UseKalman
        if(vars.currentPosition - vars.SamplesInChunk)-1 <= 0
            sample =  EEG.Recording((vars.currentPosition - vars.SamplesInChunk):vars.currentPosition - 1, EEG.PrimaryChannel);
            sample = [0; sample];
        else
            sample =  EEG.Recording((vars.currentPosition - vars.SamplesInChunk)-1:vars.currentPosition - 1, EEG.PrimaryChannel);
        end
    else
        if(vars.currentPosition - vars.SamplesInChunk)-1 <= 0
            
            sample =  EEG.Kalman_Signal((vars.currentPosition - vars.SamplesInChunk):vars.currentPosition - 1, EEG.KalmanPrimary);
            sample = [0; sample];
        else
            sample =  EEG.Kalman_Signal((vars.currentPosition - vars.SamplesInChunk)-1:vars.currentPosition - 1, EEG.KalmanPrimary);
        end
    end
    
%     [sample, vars.zhp] = filter(vars.b_hp, vars.a_hp, sample(2:end), vars.zhp);
%     sample = [0; sample];
    
    [FiltSample, vars.z] = filter(vars.b, vars.a, sample(2:end), vars.z); 
    X = zeros(3, length(sample) - 1, 1, 1);
    X(1, :, 1, 1) = sample(2:end);
    X(2, :, 1, 1) = diff(sample);
    X(3, :, 1, 1) = FiltSample;
    [vars.PhasePredictor, Pred] = predictAndUpdateState(vars.PhasePredictor, {X});
    PredAngle = angle(Pred{1}(1, end) + sqrt(-1) * Pred{1}(2, end));
%     vars.X(:, vars.currentPosition - 1) = X(:, end, 1, 1);
%     vars.Angles(vars.currentPosition - 1) = PredAngle;
 %   fprintf('hi: %f\n', PredAngle);     
    if (vars.currentPosition - vars.TriggerBuffer) > vars.LastStimPosition
        
        Mag = norm(Pred{1}(:, end));
        vars.allMags(end+1) = Mag;

        if Mag > EEG.Threshold
            if (PredAngle>=deg2rad(-60) && PredAngle<=deg2rad(-10))
                idx = (vars.currentPosition-EEG.fs*30-1):(vars.currentPosition-1);
                if idx(1) >= 1
                    
                    delp = mean(envelope(filter(vars.b_delta, vars.a_delta, EEG.Recording(idx,9)), length(idx), 'rms'));
                    vars.alldelps(end+1) = delp;

                    movs = mean(envelope(filter(vars.b_mov, vars.a_mov, EEG.Recording(idx,9)), length(idx), 'rms'));
                    vars.allmovs(end+1) = movs;
                    
                    disp('delp is ' + string(delp) + ' and movs is ' + string(movs))
                    
                    if movs < vars.movsthresh                      
                        if delp > vars.delpthresh 
                            %PsychPortAudio('Start', vars.audio_port, vars.repetitions, vars.ChunkTime + vars.SlowWaveDelay, 0);
                            %sound(Sound, fsSound)
                            vars.StimTimes(vars.StimCount) = round(vars.currentPosition + vars.SlowWaveDelay * EEG.fs);
                            vars.StimCount = vars.StimCount + 1;
                            vars.LastStimPosition = vars.currentPosition;
                            toc
                            disp('Mag is ' + string(Mag))
                            disp(vars.StimCount)
                        end
                    end
                end
            end
        end
    end
end
tic
end