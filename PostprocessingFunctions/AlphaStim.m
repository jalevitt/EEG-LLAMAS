function [vars, Graphs, EEG] = AlphaStim(EEG, vars, Graphs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if vars.UseAlphaStim && vars.currentPosition > vars.alphaPosition
    %select our sample and filter it
    sampleForAlpha = EEG.Recording(vars.currentPosition - EEG.fs + 1:vars.currentPosition, EEG.PrimaryChannel);
    FiltData = filtfilt(vars.highPassIIR, sampleForAlpha);

    %calc PSD using welch's method
    [psd, ~] = pwelch(FiltData, [], [], [], EEG.fs);

    %normalize PSD
    psd = psd/sum(psd(vars.fullMin:vars.fullMax));

    %calculate band power, and store it, and a smoothed version
    alphaIDX = vars.alphaPosition/(EEG.fs / vars.alphaPerSecond);
    vars.alphaVec(alphaIDX) = mean(psd(vars.alphaMin:vars.alphaMax));
    vars.smoothedAlphaVec(alphaIDX) = dot(vars.alphaVec(alphaIDX - 4:alphaIDX), vars.sawTooth);

    %Determine if the subject's eyes are closed
    [~, score] = predict(vars.EOECModel,  vars.smoothedAlphaVec(alphaIDX));
    vars.pEyesClosed(alphaIDX) = score(1) / sum(score);

    %if their eyes are closed, we can look for a chance to emit a sound
    if vars.UseTriggers && vars.pEyesClosed(alphaIDX) > 0.5

        %get latest quarter second of data
        FiltData = sampleForAlpha(end - EEG.fs/4:end);
        %filter it
        FiltData = filtfilt(vars.bandPassFIR, FiltData);

        %Use the mansouri phase predition algorithm to create a
        %predition 1/2 second long
        [Prediction] = MansouriPhase(FiltData', EEG.fs, EEG.fs / 2);

        % read from timer to estimate how much time has passed since
        % chunk was received
        compTime = toc(vars.clock);

        % find the delay until the next upcoming peak in our prediction, after ignoring
        % samples lost to compTime and sound prodution lag
        [~, delay] = max(Prediction(round(compTime * EEG.fs + 1 + vars.IntrinsicLag * EEG.fs):end));

        %convert delay to seconds
        delay = delay/EEG.fs;

        %Play sound with delay
        PsychPortAudio('Start', vars.audio_port, vars.repetitions, GetSecs() + delay, 0);
        t = round(vars.currentPosition + vars.IntrinsicLag * EEG.fs + compTime * EEG.fs + delay * EEG.fs);
        if ~numel(t) == 0

            vars.StimTimes(vars.StimCount) = t(1);
        else
            vars.StimTimes(vars.StimCount) = vars.currentPosition;
        end
        vars.StimCount = vars.StimCount + 1;
        vars.LastStimPosition = vars.currentPosition;

    end

    %update our graph of p(eyesClosed)
    if isgraphics(Graphs.Eyes)
        set(0, 'CurrentFigure', Graphs.Eyes)
        plot((1:length(vars.alphaVec))/vars.alphaPerSecond, vars.pEyesClosed)
        hold on
        xline(vars.currentPosition / EEG.fs, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.25);
        hold off
        xlabel('Time (s)')
        ylabel('Probability eyes closed')
        ylim([-0.5,1.5])
    end

    vars.alphaPosition = vars.alphaPosition + EEG.fs / vars.alphaPerSecond;
end


end

