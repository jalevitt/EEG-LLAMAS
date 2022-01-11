function [EEG] = OfflineGAC(EEG)
% Input is EEG struct from RTEEG 
% output is EEG struct with added OfflineGAC field.
% EEG must have the following fields:
    % EEG.trGap -  the expected samples per TR
    % EEG.vars.trMarker - the value in the event channel that corresponds
    %   with a TR
    % EEG.Recording_orig - the raw EEG data
    % EEG.numChans - the nuber of EEG channels

TRGap = EEG.trGap;
TRMark = EEG.vars.trMarker;
TRIDX = find(EEG.Recording_orig(:, EEG.numChans + 1) == TRMark);
EEG.OfflineGAC = EEG.Recording_orig;
ntr = 20;
for i = ntr + 1:length(TRIDX)
    idx = TRIDX(i);
    past = EEG.Recording_orig(idx - TRGap * ntr + 1:idx, 1:EEG.numChans);
    Template = zeros(TRGap, EEG.numChans);
    for j = 1:ntr
       Template = Template + past((j - 1) * TRGap + 1:j*TRGap, :); 
    end
    Template = Template / ntr;
    Block = EEG.Recording_orig(idx + 1 - TRGap:idx, 1:EEG.numChans);
    Block = Block - Template;
    EEG.OfflineGAC(idx + 1 - TRGap:idx, 1:EEG.numChans) = Block;
end
end

