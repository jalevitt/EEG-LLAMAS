% Remove gradient artifact with moving average.
% data: channelsxtimepoints
% trig: index of tr occurrences


function clean=removegac(data,trig,ntr);

if nargin<3
ntr=10; % use for number of TRs to average.
end
tr=max(diff(trig));
disp(sprintf('Found %d TRs with %0.2f ms separation',length(trig),tr));

clean=data;

for i=2:length(trig)
    %%
    template=zeros(size(data,1),tr);
    na=0;
    for j=max(1,i-ntr):i-1
        template=template+(data(:,trig(j):trig(j)+tr-1)')';
        na=na+1;
    end
    if trig(i)+tr<length(data)
        clean(:,trig(i):trig(i)+tr-1)=data(:,trig(i):trig(i)+tr-1)-template/na;
    else
        % find padded endpoint
        ti=trig(i):length(data);
        ti=ti(data(ti)~=0);
        clean(:,ti)=data(:,ti)-template(:,1:length(ti))/na;
    end
end




