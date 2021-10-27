% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();


chunkSize = 40;
fs =  5000;
pauseTime = chunkSize/fs;
numChans = 65;

% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Test2','EEG', numChans, fs,'cf_float32','sdfwerr32432');

disp('Opening an outlet...');
outlet = lsl_outlet(info);
disp('Beginning Transmission');
i = 1; %round(160 * fs);
%send chunks until we reach the end of the file 
while i == i
    chunk = rand(chunkSize, numChans);
    chunk(:, end) = chunk(:, end) > .9999;
    outlet.push_chunk(chunk');
    pause(pauseTime)
    i = i + chunkSize;
    
end
outlet.delete()
disp('Ending Transmission')