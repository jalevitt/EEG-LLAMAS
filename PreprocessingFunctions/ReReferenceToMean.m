function [vars, Graph, EEG] = ReReferenceToMean(EEG, vars, Graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if EEG.MeanRereference == 1
    idx = 1:EEG.numChans;
    M = mean(vars.OrigChunk(idx, :), 1);
    if ~vars.UseTriggers
        vars.OrigChunk = vars.OrigChunk - M;
    else
        vars.OrigChunk(1:end-1, :) = vars.OrigChunk(1:end-1, :) - M;      
    end
end

end

