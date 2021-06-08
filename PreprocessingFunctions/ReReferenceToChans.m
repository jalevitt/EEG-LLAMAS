function [vars, Graph, EEG] = ReReferenceToChans(EEG, vars, Graph)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if EEG.ChansRereference == 1
    if ~vars.UseTriggers   
        M = mean(vars.OrigChunk(EEG.ReReferenceChans, :), 1);
        vars.OrigChunk = vars.OrigChunk - M;
    else
        M = mean(vars.OrigChunk(EEG.ReReferenceChans, :), 1);
        vars.OrigChunk(1:end-1, :) = vars.OrigChunk(1:end-1, :) - M;
    end
end

end

