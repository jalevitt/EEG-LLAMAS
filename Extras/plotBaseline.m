figure()
plot(EEG.vars.alldelps)
title('delp')
figure()
plot(EEG.vars.allmovs)
title('movs')
figure()
plot(EEG.vars.allMags)
title('mags')

%%

params=struct;
TK = 5;
K = 2*TK-1;
params.tapers=[TK K]; % [TW K] TK= time-bandwidth = sampling frequency*spectral width, K = 2*TW - 1

params.Fs=200;

movingwin=[5 .1]; % [windowlength windowslide] 
params.fpass=[0 20]; % frequency range 

dataa = EEG.Recording';
[spec, t, f] = multitaper_spectrogram(dataa(9,:),params.Fs,params.fpass, params.tapers, movingwin);

%%
medmag = median(EEG.vars.allMags)

temp=sort(EEG.vars.allMags,'descend');
top5p = median(temp(1:floor(length(temp)/20)))
top2p = median(temp(1:floor(length(temp)/50)))
top1p = median(temp(1:floor(length(temp)/100)))

%%

poolobj = gcp('nocreate');
delete(poolobj);
