% Load the data
cd('raw/run01')
list = ls();
load('hdr.mat')
load('ecg.mat')
ecg = ch;
load('emg.mat')
emg = ch;
load('resp.mat')
resp = ch;
eeg = zeros(100, length(ecg));
for f = 1:size(list, 1)
    if strcmp(list(f, 1), 'c')
        splt = split(list(f, :), '.');
        idx = str2num(splt{1}(2:end));
        load(list( f, :))
        eeg(idx, :) = ch;
    end
end
cd('../..')

%%

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
oldfs = orighdr.Fs;
newfs = 200;
ds = oldfs / newfs;

filtsize = 2000;
ff = firls(filtsize, [0 50 55 oldfs / 2] / oldfs * 2, [1 1 0 0]);
x = filter(ff, 1, data, [], 2);
x = circshift(x,[0 -filtsize / 2]);

x = x(:, floor(ds / 2):ds:end);


%% gradient artifact cleaning
gradtr=20;
cleandata = x;
trsamp_ds = floor(hdr.trsamp / ds);
if length(unique(diff(trsamp_ds))) > 1
    warning('WARNING: TR times are uneven');    
end

cleandata = removegac(cleandata, trsamp_ds, gradtr);
cleandata_orig = removegac(data, trsamp, gradtr);
x_orig = filter(ff, 1, cleandata_orig, [], 2);
x_orig = circshift(x_orig,[0 -filtsize / 2]);

x_orig = x_orig(:, floor(ds / 2):ds:end);
%%
n = size(x, 2);
n_old = size(eeg, 2);
ts = linspace(1/newfs, n / newfs, n);
oldts = linspace(1/oldfs, n_old / oldfs, n_old);
figure
plot(ts, x(10, :), ts, cleandata(10, :), ts, x_orig(10, :))
xlabel('time(s)')
legend({'before', 'after', 'orig'})
hold on
% for i=1:length(trsamp)
%     xline(trsamp(i) / newfs);
% end