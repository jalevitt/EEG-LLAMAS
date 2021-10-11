function [EEG] = LoadMFF(fname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath('C:\Users\lewislab\Documents\Matlab\fieldtrip-20200821');
ft_defaults;

data=ft_read_data(fname);
orighdr=ft_read_header(fname);
event=ft_read_event(fname);

evflag = 1;
hdr=struct;
hdr.chantype=orighdr.chantype([1:256 258:end]);
nc=length(hdr.chantype);
np=nc-256;
hdr.refchan=cell(nc,1);
hdr.chantype(1:224)={'eeg'};
hdr.chantype(225:256)={'cheek'};
hdr.refchan(1:256)={'cz'};
hdr.refchan(257:end)={'bipolar'};
hdr.processlog={};
hdr.fband=[];
hdr.badchans=[];
hdr.events=event; % store event structure
hdr.fs=orighdr.Fs;
if evflag&any(strcmp({event.value},'TREV'))
    trcase=1;
    allev= round([event.sample]);
    trs=allev(strcmp({event.value},'TREV'));
    if length(unique(diff(trs))>1)
        di=diff(trs);
        m=mode(di);
        starts=find([1 diff(di) 1]);
        [m mi]=max(diff(starts));
        tt=trs(mi:starts(mi+1));
        trs=tt;
    end    
    hdr.trsamp=trs; % tr in samples
    hdr.trs=hdr.trsamp/hdr.fs; % tr in seconds
    trsamp=hdr.trsamp;
else 
    trcase=0;
end
if evflag
hdr.evtime=[hdr.events.sample]/hdr.fs;
hdr.evlabel={hdr.events.value};
else
    hdr.evtime=[];
    hdr.evlabel={};
end
labels=cell(256,1);
for i=1:length(labels)
    labels{i}=['c' num2str(i)];
end
hdr.label=labels;


fs = hdr.fs;
data = double(data);
l = size(data, 2);
Triggers = zeros(1, l);
Triggers(trsamp) = 1;
data = [data; Triggers];

EEG.Recording = data';
EEG.fs = fs;
EEG.hdr = hdr;

end

