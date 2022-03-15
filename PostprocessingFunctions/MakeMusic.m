function [vars, Graph, EEG] = MakeMusic(EEG, vars, Graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(vars, 'MusicVar')
    vars.MusicVar.MusicFlag = 1;

    vars.MusicVar.MusicChunkSize = round(EEG.fs * 0.5);
    vars.MusicVar.MusicWindowSize = round(EEG.fs * 0.75);
    vars.MusicVar.LastMusicChunk = 200;
    vars.MusicVar.FreqData = zeros(100000, 2);
    vars.MusicVar.MusicFs = 44100;
    vars.MusicVar.NoteDur = round(0.5 * vars.MusicVar.MusicFs);
    x = (1:vars.MusicVar.NoteDur)/vars.MusicVar.MusicFs;
    vars.MusicVar.Sounds{1} = sin(x * 2 * pi * 146.83) + sin(x * 2 * pi * 196);
    vars.MusicVar.Sounds{2} = sin(x * 2 * pi * 174.61) + sin(x * 2 * pi * 246.94);
    vars.MusicVar.Sounds{3} = sin(x * 2 * pi * 196) + sin(x * 2 * pi * 261.63);
    vars.MusicVar.Sounds{4} = sin(x * 2 * pi * 207.65) + sin(x * 2 * pi * 277.18);
    vars.MusicVar.AudioWriter = audioDeviceWriter('SampleRate', vars.MusicVar.MusicFs);
    vars.MusicVar.MusicCount = 1;
    vars.MusicVar.Freqs = [6, 10, 15, 20];
elseif vars.currentPosition > vars.MusicVar.LastMusicChunk + vars.MusicVar.MusicChunkSize
    fxx = 0:22;
    musicChunk = EEG.Recording((vars.currentPosition - vars.MusicVar.MusicWindowSize):vars.currentPosition - 1 , EEG.PrimaryChannel);
    chunkSpec = pwelch(musicChunk, 100, 25, fxx, EEG.fs);
    chunkSpec = mag2db(chunkSpec);
    
    b = robustfit(fxx, chunkSpec);
    Residual = chunkSpec - (fxx) * b(2) + b(1);
    if isgraphics(Graph.Eyes)
        set(0, 'CurrentFigure', Graph.Eyes)
        plot(fxx, chunkSpec, fxx, Residual)
    end
    [~, idx]  = max(Residual);
    freq = fxx(idx);
    [~, idx] = min(abs(vars.MusicVar.Freqs - freq));
    closest = vars.MusicVar.Freqs(idx);
    if closest == vars.MusicVar.Freqs(1)
        toPlay = vars.MusicVar.Sounds{1};
    elseif closest == vars.MusicVar.Freqs(2)
         toPlay = vars.MusicVar.Sounds{2};
    elseif closest == vars.MusicVar.Freqs(3)
         toPlay = vars.MusicVar.Sounds{3};
    elseif closest == vars.MusicVar.Freqs(4)
         toPlay = vars.MusicVar.Sounds{4};
    end
    PsychPortAudio('FillBuffer', vars.audio_port, [toPlay; toPlay]);
    PsychPortAudio('Start', vars.audio_port, 1, 0);
%     for j = 1:1024:length(toPlay)
%        if j + 1024 > length(toPlay)
%            sub = [toPlay(j:end); zeros(1024 - (1 + length(toPlay) - j), 1)];
%        else
%            sub = toPlay(j:j + 1024 - 1);
%        end
%        vars.MusicVar.AudioWriter(sub');
%     end
    vars.MusicVar.FreqData(vars.MusicVar.MusicCount, 1) = vars.currentPosition;
    vars.MusicVar.FreqData(vars.MusicVar.MusicCount, 2) = freq;
    vars.MusicVar.MusicCount = vars.MusicVar.MusicCount + 1;
    fprintf('The peak frequency is %d \n', freq)
    vars.MusicVar.LastMusicChunk = vars.MusicVar.LastMusicChunk + vars.MusicVar.MusicChunkSize;
    
end




end

