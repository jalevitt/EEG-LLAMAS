function [responseRate, meanERPs] = AnalyzeSoundTest(EEG)
% analyze the results of a SoundTest run
%   Detailed explanation goes here

if ~EEG.vars.UseKalman
    data = EEG.Recording(:, EEG.PrimaryChannel);
else
    data = EEG.Kalman_Signal(:, EEG.KalmanPrimary);
end

numStims = EEG.vars.SoundTest.StimCount;
stimTimes = EEG.vars.SoundTest.StimTimes(1:numStims);
stimVols = EEG.vars.SoundTest.StimVols(1:numStims);
trig = EEG.vars.BtnMarker;
trigTimes = zeros(size(stimTimes));
allTrigs = find(EEG.Recording(:, end - 1) == trig);
behaveResponse = zeros(size(stimTimes));



% get subject behavior
[row, column, dataValues] = find(EEG.KeyPresses);



responseCols = [11, 12];
trCols = [15];

firstTrig = find(EEG.Recording(:, end - 1) == EEG.vars.trMarker, 1, 'first');
trTimes = dataValues(column == trCols);
buttonTimes = (dataValues - min(trTimes)) * EEG.fs + firstTrig;

verifiedTimes = [];
for i = 1:length(responseCols)
    verifiedTimes = [verifiedTimes; buttonTimes(column == responseCols(i))];
end
verifiedTimes = sort(verifiedTimes);


% find the trigger aand behav corresponding to each stim

ERPstart = -50;
ERPend = 150;
ERPs = zeros(numStims - 1, 1 + ERPend - ERPstart);
for i = 1:numStims - 1
    time = stimTimes(i);
    nextTime = stimTimes(i + 1);
    f = find((allTrigs > time) & (allTrigs < nextTime));
    if isempty(f)
        trigTimes(i) = nan;
    else
        trigTimes(i) = allTrigs(min(f));
    end  
    f = find((verifiedTimes > time) & (verifiedTimes < nextTime));
    if ~isempty(f)
        behaveResponse(i) = 1;
    end
    ERPs(i, :) = data(time + ERPstart:time + ERPend); 
end


responseRate = zeros(size(EEG.vars.SoundTest.VolOptions));
meanERPs = zeros(length(responseRate), 1 + ERPend - ERPstart);
for i = 1:length(responseRate)
   num_this_vol = sum(stimVols(1:end - 1) == EEG.vars.SoundTest.VolOptions(i));
   idx_this_vol = find(stimVols(1:end - 1) == EEG.vars.SoundTest.VolOptions(i));
   numResponse = sum(behaveResponse(idx_this_vol));
   responseRate(i) = 100 * numResponse / num_this_vol;  
   meanERP(i, :) = mean(ERPs(idx_this_vol,:), 1);
end

t = (ERPstart:ERPend) / EEG.fs;

figure
bar(responseRate)
xlabel('Volume')
xticks(1:length(responseRate))
xticklabels(EEG.vars.SoundTest.VolOptions)
ylabel('Response Rate (%)')
ylim([0, 100])
xlim([0, length(responseRate) + 1])

figure
for i = 1:length(responseRate)
    subplot(2, 5, i)
    
    plot(t, meanERP(i, :))
    xlabel('time (s)')
    ylabel('amplitude (uV)')
    ylim([-20, 20])
    xlim([ERPstart, ERPend]/EEG.fs)
    title(sprintf('Volume = %f', EEG.vars.SoundTest.VolOptions(i)))
end



end

