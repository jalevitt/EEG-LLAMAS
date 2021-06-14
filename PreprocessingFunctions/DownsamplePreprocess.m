function [vars, Graphs, EEG] = DownsamplePreprocess(EEG, vars, Graphs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if vars.UseTriggers
   [O, vars.zi] = filter(vars.DS_Lowpass, 1, vars.OrigChunk(1:end - 1, :)', vars.zi, 1);
   vars.OrigChunk(1:end - 1, :) = O';
else
    [O, vars.zi] = filter(vars.DS_Lowpass, 1, vars.OrigChunk', vars.zi, 1);
    vars.OrigChunk = O';
end
%down sampling works slightly differently if we're using
%triggers
if vars.UseTriggersvars
    vars.stamps = downsample(vars.stamps', vars.DSrate, vars.dsBuffer)';

    %triggers must be downsampled using this specialized
    %function
    [vars.chunk, vars.dsBuffer] = DownSampleTriggrs(vars.OrigChunk, vars.DSrate, vars.dsBuffer);

    [vars.SamplesInChunk, vars.ChansInChunk] = size(vars.chunk');
    if vars.SamplesInChunk ~= length(vars.stamps)
        vars.stamps = zeros(1, vars.SamplesInChunk);
    end
else

    %data without triggers can be downsampled normally
    if vars.SamplesInChunk == 1 && (vars.dsBuffer == 0)
        vars.chunk = vars.OrigChunk;
    elseif vars.SamplesInChunk == 1
        vars.chunk = zeros(vars.ChansInChunk, 0);
    else
        vars.chunk = downsample(vars.OrigChunk', vars.DSrate, dsBuffer)';
    end
    vars.stamps = downsample(vars.stamps', vars.DSrate, vars.dsBuffer)';
    vars.dsBuffer = mod(vars.SamplesInChunk - vars.dsBuffer, vars.DSrate);
    [vars.SamplesInChunk, ~] = size(vars.chunk');
end
end

