addpath('C:\Users\lewislab\Documents\Matlab\fieldtrip-20200821');
ft_defaults;

loadfile = dir('C:\Users\lewislab\Desktop\RealTimeEEGTestData\osceeg\mff');

fname = 'C:\Users\lewislab\Desktop\RealTimeEEGTestData\osceeg\mff\realrun1_20161101_100049.mff';
data=ft_read_data(fname);
orighdr=ft_read_header(fname);
event=ft_read_event(fname);
%%
bcg=1:224;
border=[91 82 73 67 61 220 219 218 217 216];
eeg=[79 78 101 126 150 137 136 116 96 95 72 52 41 21 15 62 38 143 183 207 24 221 204 193 172 160 168 37 32 25 18 10 31];
bcg=bcg(~ismember(bcg,eeg));
border=border(~ismember(border,eeg));
%% Set variables
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

%%
chunkSize = 32;
fs = hdr.fs;
FudgeFactor = 25/31;
pauseTime = FudgeFactor * chunkSize/fs;
data = double(data);
[numChans, l] = size(data);
Triggers = zeros(1, l);
Triggers(trsamp) = 1;
data = [data; Triggers];
% make a new stream outlet
disp('Loading the library...');
lib = lsl_loadlib();

disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Test2','EEG', numChans + 1, fs,'cf_float32','sdfwerr32432');

disp('Opening an outlet...');
outlet = lsl_outlet(info);
disp('Beginning Transmission');
i = 1; %round(160 * fs);
%send chunks until we reach the end of the file 
while i < l - chunkSize
    outlet.push_chunk(data(:, i:i+chunkSize - 1));
    pause(pauseTime)
    i = i + chunkSize
    
end
outlet.delete()
disp('Ending Transmission')
