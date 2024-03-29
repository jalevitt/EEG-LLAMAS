function [] = WriteLogInfo(handles, fid, txtdate, vars, EEG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fprintf(fid, 'Recording Name: %s\n', EEG.RecordingName');
fprintf(fid, 'Start Time: %s\n', txtdate);
fprintf(fid, 'Native Sample Rate: %d Hz\n', str2num(handles.NativeFS.String));
fprintf(fid, 'Number of channels: %d\n', str2num(handles.numChans.String));
if handles.UseDownsampling.Value
    fprintf(fid, 'Downsampling: True\n');
    fprintf(fid, 'Downsampled rate: %d\n', str2num(handles.TargetFS.String));
else
    fprintf(fid, 'Downsampling: False\n');
end
if vars.ReReferenceToChans
    fprintf(fid, 'Re-Referenced to Channels: '); 
    for i = 1:length(EEG.ReReferenceChans)
        fprintf(fid, '%d, ', EEG.ReReferenceChans(i)); 
    end
    fprintf(fid, '\n'); 
else
    fprintf(fid, 'Native referencing Montage\n ');
end
if handles.UseTriggers.Value
    fprintf(fid, 'Triggers: True\n');
else
    fprintf(fid, 'Triggers: False\n');
end
fprintf(fid, 'Primary Channel: %d\n', str2num(handles.PrimaryChannel.String));
if vars.UseAlphaStim
    fprintf(fid, 'Alpha Stimulation: True\n');
else
    fprintf(fid, 'Alpha Stimulation: False\n');
end
if vars.UseSlowWaveStim
    fprintf(fid, 'Slow Wave Stimulation: True\n');
    fprintf(fid, 'Stimulation Threshold: %0.2f\n', str2num(handles.StimThresh.String));
else
    fprintf(fid, 'Slow Wave Stimulation: False\n');
end
if handles.GAC.Value
    fprintf(fid, 'Gradient Artifact Correction: True\n');
    fprintf(fid, '\tntr = %d \n', vars.ntr);
else
    fprintf(fid, 'Gradient Artifact Correction: False\n');
end
if handles.UseOld.Value
    fprintf(fid, 'Is Replay: True\n');
    fprintf(fid, 'Path to Replayed data: %s\n', handles.ReplayPath);
else
    fprintf(fid, 'Is Replay: False\n');
end
if vars.UseKalman
    fprintf(fid, 'Use Kalman Filter: True\n');
else
    fprintf(fid, 'Use Kalman Filter: False\n');
end    
end

