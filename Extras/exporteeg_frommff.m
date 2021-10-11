% exporteeg_frommff:

% structure to be exported:
% file for each channel, containing 'data','fs','label'
% file for whole run, containing:
% 'trs','events','ts','badchans','fband','refchan',
% 'processlog','runname','chlabels','chantype'

% this creates raw file and filtered file; later, go mark bad channels both
% manually and automatically

function [x,hdr]=exporteeg_frommff(subj,runname,ploton)
if nargin<3
    ploton=0;
end

disp(['starting ' subj runname]);
pathstem='/cluster/ldl/ldlewis/eegfmri/';
loadfile=dir([pathstem subj '/mff/*' runname '*']);
if isempty(loadfile)
    loadfile=dir([pathstem subj '/mff/*' runname(1:end-2) runname(end) '*']);
end
savepath=[pathstem subj '/eeg/raw/' runname '/'];
mkdir(savepath);

fname=[pathstem subj '/mff/' loadfile(1).name];

disp(['reading file: ' fname]);
orighdr=ft_read_header(fname);
data=ft_read_data(fname);

try
event=ft_read_event(fname);
evflag=1;
catch
    disp('no events found');
    event=[];
    evflag=0;
end

%% Set variables

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
    allev=[event.sample];
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

% label physio channels
switch np
    case 1
        hdr.chantype{end}='ecg';
        hdr.label{end+1}='ecg';
    case 3
         hdr.chantype{end-2}='emg';
        hdr.chantype{end-1}='ecg';
        hdr.chantype{end}='resp';
        hdr.label{end+1}='emg';
        hdr.label{end+1}='ecg';
        hdr.label{end+1}='resp';
end
disp('warning: check labels for physio, set to:');
hdr.label(end-np+1:end)

% manually enter bcg channels!!
[bcgchans,borderchans]=getbcgchans(subj);
hdr.chantype([bcgchans])={'bcg'};
hdr.chantype([borderchans])={'border'};

% key runwise vars
ts=(1:length(data))/hdr.fs;
hdr.ts=ts;

%% gradient artifact cleaning
gradtr=20;
cleandata=data([1:256 258:end],:);
if trcase
    disp(['Cleaning gradient artifacts with n=' num2str(gradtr)]);
        if length(unique(diff(trsamp)))>1
        error('WARNING: TR times are uneven');
        1      
    end
    % clean all channels except respiratory % changed march 30, 201: also
    % clean respiratory! but with more TRs to avoid any highpass filtering
    cc=(~strcmp(hdr.chantype,'resp'));
    respchan=find((~strcmp(hdr.chantype,'resp'))); %  added march 30 2017
    cleandata(cc,:)=removegac(cleandata(cc,:),trsamp,gradtr);
    cleandata(respchan,:)=removegac(cleandata(respchan,:),trsamp,40);%  added march 30 2017
    hdr.processlog{end+1}=['remTR_' num2str(gradtr)];
    if ploton
        figure;
        plot(diff(trsamp));
        title('TR difference');
    end

end


%% save files
if size(cleandata,1)~=length(hdr.label)
    error('check size')
end
disp(['Saving into ' savepath ]);
save([savepath 'hdr'],'hdr');
for c=1:size(cleandata,1)
    ch=cleandata(c,:);
    save([savepath hdr.label{c}],'ch');
end

%% Check if there are impedances stored.
calibtypes={};
for i=1:length(orighdr.orig.xml.info1.calibrations)
    calibtypes{i}=orighdr.orig.xml.info1.calibrations(i).calibration.type;
end
if any(strcmp(calibtypes,'ICAL'))
    ix=find(strcmp(calibtypes,'ICAL'));
    ix=ix(end);
    impedances=zeros(257,1);
    for j=1:257
        impedances(j)=str2num(orighdr.orig.xml.info1.calibrations(ix).calibration.channels(j).ch.ch);
    end
    disp(['saving impedances as ' savepath '../../impedances_' runname '.mat']);
    save([savepath '../../impedances_' runname '.mat'],'impedances');
end

%% filter and downsample for most uses
oldfs=hdr.fs;
newfs=200;

filtsize=2000;
ff=firls(filtsize,[0 50 55 oldfs/2]/oldfs*2,[1 1 0 0]);
emgff=firls(filtsize,[0 70 75 oldfs/2]/oldfs*2,[1 1 0 0]);
ds=oldfs/newfs;

% filter all channels
x=filter(ff,1,cleandata,[],2);
x=circshift(x,[0 -filtsize/2]);
% if there is an EMG channel, filter it differently
emg=find(strcmp(hdr.chantype,'emg'));
if ~isempty(emg)
    x(emg,:)=filter(emgff,1,cleandata(emg,:),[],2);
    x(emg,:)=circshift(x(emg,:),[0 -filtsize/2]);
end

x=x(:,floor(ds/2):ds:end);
oldts=hdr.ts;
newts=hdr.ts(floor(ds/2):ds:end);
hdr.ts=newts;
hdr.fs=newfs;
hdr.processlog{end+1}='filter EEG 0-50hz (EMG 0-70 hz)';
hdr.fband=[0 50];

savepath=[pathstem subj '/eeg/fs200/' runname '/'];
mkdir(savepath);
disp(['Saving filtered data into ' savepath ]);
save([savepath 'hdr'],'hdr');
for c=1:size(x,1)
    ch=x(c,:);
    save([savepath hdr.label{c}],'ch');
end

%% plot to make sure everything looks reasonable
if ploton
    figure
    plot(oldts,cleandata(126,:),newts,x(126,:))%,hdr.ts,cleandata(257,:));
    legend({'orig','filt'});
    xlabel('Time (s)');
end




