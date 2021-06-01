
%build the sound we'll play
fsSound = vars.audioDevices(vars.deviceID + 1).DefaultSampleRate; %48000;
s = 0.05;
cn = dsp.ColoredNoise(1, fsSound * s, 1, 'OutputDataType', 'double');
PinkSound = cn();
PinkSound = hamming(length(PinkSound)).*PinkSound;

%initialize the audio player
InitializePsychSound;

vars.audioDevices = PsychPortAudio('GetDevices');
if isunix
    vars.deviceID = 0;
else
    vars.deviceID = 5;
end

audio_port = PsychPortAudio('Open', vars.deviceID, 1, 3, fsSound, 2, []);
audio_to_play = [PinkSound, PinkSound]';
PsychPortAudio('FillBuffer', audio_port, audio_to_play);

% Play the Sound
PsychPortAudio('Start', audio_port, 1, GetSecs(), 0);

pause(0.5)
%clean up
PsychPortAudio('Close')