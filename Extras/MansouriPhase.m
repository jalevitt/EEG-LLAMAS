function [Prediction] = MansouriPhase(signal, fs, nPoints)
%Takes ashort EEG signal, sample rate, and nPoints as input,
%Returns a projected signal nPoints in length based on Mansouri et al 2017
signal = signal - mean(signal);
L = length(signal);
L_adjusted = 10000;
FreqDomain = fft([signal, zeros(1, L_adjusted - L)]);
[~, mxBin] = max(abs(FreqDomain(1:L_adjusted/2)));

phase = angle(FreqDomain(mxBin));
freq = (mxBin - 1) * fs/L_adjusted;

time = linspace(0, (L + nPoints - 1)/fs, nPoints + L);
time(1:L) = [];
Prediction = cos(time * 2 * pi * freq + phase);

end

