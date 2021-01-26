function [vars] = SlowWaveStim(EEG, vars)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if vars.UseSlowWaveStim && vars.SamplesInChunk > 0 
    if ~vars.UseKalman
        sample =  EEG.Recording(vars.currentPosition - vars.SamplesInChunk:vars.currentPosition - 1, EEG.PrimaryChannel);
    else
        sample =  EEG.Kalman_Signal(vars.currentPosition - vars.SamplesInChunk:vars.currentPosition - 1, EEG.PrimaryChannel);
    end
    [FiltSample, vars.z] = filter(vars.b, vars.a, sample, vars.z);

    if max(FiltSample) > EEG.Threshold  && (vars.currentPosition - vars.TriggerBuffer) > vars.LastStimPosition
        PsychPortAudio('Start', vars.audio_port, vars.repetitions, GetSecs() + vars.SlowWaveDelay, 0);
        %sound(Sound, fsSound)
        vars.StimTimes(vars.StimCount) = vars.currentPosition;
        vars.StimCount = vars.StimCount + 1;
        vars.LastStimPosition = vars.currentPosition;
    end
end
end

