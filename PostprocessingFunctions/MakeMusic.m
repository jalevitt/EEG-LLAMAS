function [vars, Graph, EEG] = MakeMusic(EEG, vars, Graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(vars, 'MusicFlag')
    vars.MusicFlag = 1;

    vars.MusicChunkSize = round(EEG.fs * 0.25);
    vars.MusicWindowSize = round(EEG.fs * 1);
    vars.LastMusicChunk = 200;
elseif vars.currentPosition > vars.LastMusicChunk + vars.MusicChunkSize
    musicChunk = EEG.Recording((vars.currentPosition - vars.MusicWindowSize):vars.currentPosition - 1 , EEG.PrimaryChannel);
    chunkSpec = pwelch(musicChunk, 100, 25, 1:25, EEG.fs);
    chunkSpec = mag2db(chunkSpec);
    
    b = robustfit(1:25, chunkSpec);
    Residual = chunkSpec - (1:25) * b(2) + b(1);
    fxx = 1:25;
    if isgraphics(Graph.Eyes)
        set(0, 'CurrentFigure', Graph.Eyes)
        plot(fxx, chunkSpec, fxx, Residual)
    end
    [~, freq]  = max(Residual);
    fprintf('The peak frequency is %d \n', freq)
    vars.LastMusicChunk = vars.LastMusicChunk + vars.MusicChunkSize;
    
end




end

