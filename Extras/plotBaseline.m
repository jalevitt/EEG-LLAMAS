figure()
plot(EEG.vars.alldelps)
figure()
plot(EEG.vars.allmovs)
figure()
plot(EEG.vars.allMags)

medmag = median(EEG.vars.allMags)

temp=sort(EEG.vars.allMags,'descend');
top5p = median(temp(1:floor(length(temp)/20)))
top2p = median(temp(1:floor(length(temp)/50)))
top1p = median(temp(1:floor(length(temp)/100)))