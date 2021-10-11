function [vars, Graph, EEG] = ReReferenceToMean(EEG, vars, Graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if EEG.MeanRereference == 1
    idx = 1:EEG.numChans;
    M = mean(vars.OrigChunk(idx, :), 1);
    vars.OrigChunk = vars.OrigChunk - M;
end

end

