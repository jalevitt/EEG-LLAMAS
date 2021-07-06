% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();



%Read in the test data
uiopen("*.mat");

chunkSize = 32;
fs = EEG.fs;
pauseTime = chunkSize/fs;
data = double(EEG.Recording(:, 1:end - 2)');
[numChans, l] = size(data);

% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Test2','EEG', numChans, fs,'cf_float32','sdfwerr32432');

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
