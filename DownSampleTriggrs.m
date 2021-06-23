function [chunk, dsBuffer] = DownSampleTriggrs(chunk, DSrate, dsBuffer)
% Downsamples a chunk of a trigger stream, while preserving all the
% positive triggers
[SamplesInChunk, ChansInChunk] = size(chunk');
triggers = chunk(ChansInChunk, :) > 0;
if sum(triggers > 0)
    idx = find(triggers);
    trigVals = chunk(ChansInChunk, idx);
    idx = idx - dsBuffer;
    idx = round(idx/DSrate);
    
    idx(idx < 1) = 1;
end

if SamplesInChunk == 1 && dsBuffer == 0 %% (dsBuffer == 0 || sum(triggers) > 0)
    chunk = chunk;
elseif SamplesInChunk == 1
    chunk = zeros(ChansInChunk, 0);
else
    chunk = downsample(chunk', DSrate, dsBuffer)';
end
dsBuffer = mod(SamplesInChunk - dsBuffer, DSrate);
[SamplesInChunk, ChansInChunk] = size(chunk');
% if rand < 0.01
%     SamplesInChunk
% end
if sum(triggers) > 0 
    chunk(ChansInChunk, :) = zeros(1, SamplesInChunk);
    chunk(ChansInChunk, idx) = trigVals;  
end

end
