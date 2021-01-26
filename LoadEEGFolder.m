function [EEG] = LoadEEGFolder(Path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cd(Path);
files = ls();
EEG.data = zeros(1);
[a, b] = size(EEG.data);
ChanID = zeros(300, 1);
for f = 1:length(files)
    disp(files(f, :))
    if strcmp(files(f, 1), 'c')
       Channel  =  load(files(f, :));
       strs = split(files(f, :), '.');
       ChannelNum = str2num(strs{1}(2:end));
       if length(Channel.ch) > b
           EEG.data = [EEG.data, zeros(a, length(Channel.ch) - b)];
           b = length(Channel.ch);
       end
       if ChannelNum > a
           EEG.data = [EEG.data; zeros(ChannelNum - a, b)];
           a = ChannelNum;
       end
       EEG.data(ChannelNum, :) = Channel.ch;
       ChanID(ChannelNum) = ChannelNum; 
    elseif strcmp(files(f, 1:3), 'hdr')
        S = load(files(f, :));
        EEG.hdr = S.hdr;
        EEG.srate = S.hdr.fs;
    end
end

% DeletionList = zeros(1, a) == ones(1, a);
% 
% for i = 1:a
%    if ~strcmp(EEG.hdr.chantype{i}, 'eeg')
%        DeletionList(i) = 1;
%    end
% end
% EEG.data(DeletionList, :) = [];
% EEG.ChanID = ChanID(~DeletionList);
[numchan, ~] = size(EEG.data);
EEG.nbchan = numchan;
cd('C:/Users/lewislab/Documents/labstreaminglayer/LSL/liblsl-Matlab')
end

