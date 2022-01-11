function [PSD, fx, t, ax] = Spec_detrend(Signal, window, overlap, fxx, fs, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if length(window) == 1
    l_win = window;
else
    l_win = length(window);
end

l_Sig = length(Signal);
[psd, fx] = periodogram(Signal(1:1 + l_win - 1), [], fxx, fs);
t = (1:(l_win - overlap):(l_Sig-l_win))/fs;
n_win = length(t);
PSD = zeros(length(psd), n_win);


i = 1;
count = 1;
while i < (l_Sig - l_win) + 1
    Win = Signal(i:i + l_win - 1);
    if length(window) > 1
        Win = Win.*window;
    end
    if nargin == 5
        Win = detrend(Win);
    else
        Win = detrend(Win, varargin{1});
    end
    [psd, fx] = periodogram(Win, [], fxx, fs);
    psd(psd == 0) = min(psd(psd ~= 0)) * .01;
    PSD(:, count) = psd;

    count = count + 1;
    i = i + l_win - overlap;
end


[X, Y] = meshgrid(t, fx);
PSD = mag2db(PSD);
ax = pcolor(X, Y, PSD);
ax.EdgeColor = 'none';
m = nanmean(PSD(:));
sd = nanstd(PSD(:));
colormap jet
caxis([m - 2 * sd, m + 2 * sd]);
h = colorbar;
set(get(h,'label'),'string','Power (dB)');
xlabel('Time (s)')
ylabel('Frequency (Hz)')
end

