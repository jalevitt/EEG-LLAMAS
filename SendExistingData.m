% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();



%Read in the test data
%load('C:/Users/lewislab/Desktop/RealTimeEEGTestData/JoshEEG.mat')
%EEG = LoadEEGFolder('C:/Users/lewislab/Desktop/RealTimeEEGTestData/base01_8b');
%load('C:\Users\lewislab\Desktop\RealTimeEEGTestData\fNIRS_Pilot\081219\run03')
%EEG = LoadRealTimeFile('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\TestEEG1_10_2_19_jeff.mat');
EEG = LoadRealTimeFile('C:\Users\lewislab\Desktop\RealTimeEEGTestData\Pilot\Test_11_20_19_Sleep1.mat');

chunkSize = 32;
fs = EEG.srate
FudgeFactor = 25/31;
pauseTime = FudgeFactor * chunkSize/fs;
data = double(EEG.data(1:end - 2, :));
[numChans, l] = size(data)

% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Test2','EEG', EEG.nbchan, fs,'cf_float32','sdfwerr32432');

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
