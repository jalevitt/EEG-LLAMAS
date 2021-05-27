function [vars, Graph, EEG] = ReReferenceToChans(EEG, vars, Graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if EEG.ChansRereference == 1
    M = mean(vars.OrigChunk(EEG.ReReferenceChans, :), 1);
    vars.OrigChunk = vars.OrigChunk - M;
end

end

