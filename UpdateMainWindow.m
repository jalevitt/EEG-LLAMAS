function [vars, Graph] = UpdateMainWindow(EEG, vars, Graph)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%determine if its time to update our plot
if vars.currentPosition > vars.EEGPlotPosition
    %print chunk size once in while so you can monitor it in real time
    %Ideally shoud be in the single digits or low teens but
    %spikes up to ~25 aren't unusual

    vars.SamplesInChunk
    %calc limits of graph

    xMin = vars.WindowSize * floor((vars.currentPosition - 1)/(vars.WindowSize*EEG.fs));
    xMax = xMin + vars.WindowSize;

    %choose samples to plot
    SampleMin = xMin * EEG.fs + 1;
    SampleMax = vars.currentPosition - 1;

    % make our time vector
    time = linspace(SampleMin/EEG.fs, SampleMax/EEG.fs, SampleMax - SampleMin + 1);


    %choose the sample section we'll be plotting
    if vars.UseKalman
        SampleToPlot = EEG.Kalman_Signal(SampleMin:SampleMax, :);
        if vars.UseTriggers
            SampleToPlot = [SampleToPlot, EEG.Recording(SampleMin:SampleMax, end - 1)];
        end
    else
        SampleToPlot = EEG.Recording(SampleMin:SampleMax, vars.ChannelsToPlot');
    end

    %adjust our sample to create vertical channel offsets
    for i = 1:vars.numChannelsToPlot
        mx = max(SampleToPlot(:, i));
        mn = min(SampleToPlot(:, i));
        SampleToPlot(:, i) = 2 * (SampleToPlot(:, i) - mn)/(mx - mn) + i * 2;
    end


    %make our graph
    if isgraphics(Graph)
        set(0, 'CurrentFigure', Graph)
        plot(time, SampleToPlot, 'k' )
        xLines = unique(round(xMin:xMax));
        hold on
        for xl = xLines
            xline(xl, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.25);
        end
        hold off
        xlim([xMin, xMax])
        ylim([1, vars.numChannelsToPlot * 2 + 2])
        yticks(1 + (2:2:vars.numChannelsToPlot * 2));
        if ~vars.UseKalman
            yticklabels(vars.ChannelNames(vars.ChannelsToPlot))
        end
        xlabel('Time (S)')
    end

    %set the next time we'll update our plot
    vars.EEGPlotPosition = vars.EEGPlotPosition + EEG.fs/vars.EEGPlotPerSecond;
end
end

