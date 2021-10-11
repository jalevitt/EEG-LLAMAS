function [vars, Graph, EEG] = UpdateMainWindow(EEG, vars, Graph)
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
    SampleToPlot = EEG.Recording(SampleMin:SampleMax, vars.ChannelsToPlot');
    
    if vars.UseKalman
        for i = 1:vars.numChannelsToPlot
            temp = ismember(vars.KalmanTargets, vars.ChannelsToPlotdIDX(i));
            if sum(temp)
                
                SampleToPlot(:, i) = EEG.Kalman_Signal(SampleMin:SampleMax, temp);
            end
        end
%         if vars.UseTriggers
%             SampleToPlot = [SampleToPlot, EEG.Recording(SampleMin:SampleMax, end - 1)];
%         end
%     else
%         
    end
    
    ChanStr = vars.ChannelNames(vars.ChannelsToPlot);
    YTickPos  = zeros(vars.numChannelsToPlot * 3, 1);
    YTickLab  = cell(vars.numChannelsToPlot * 3, 1);
    %adjust our sample to create vertical channel offsets
    for i = 1:vars.numChannelsToPlot
        if vars.ScaleSize == 0 
            mx = max(SampleToPlot(:, i));
            mn = min(SampleToPlot(:, i));
            SampleToPlot(:, i) = 2 * (SampleToPlot(:, i) - mn)/(mx - mn) + i * 2;
            YTickPos(i * 3 - 1) = i * 2 + 1;
            YTickPos(i * 3 - 2) = i * 2 + 0.5;
            YTickPos(i * 3) = i * 2 + 1.5;
            YTickLab(i * 3 - 1) = ChanStr(i);
            YTickLab(i * 3 - 2) = {sprintf('%.2f', 0.75 * mn + 0.25 * mx)};
            YTickLab(i * 3) = {sprintf('%.2f', 0.25 * mn + 0.75 * mx)};
        else
            SampleToPlot(:, i) = 2 * SampleToPlot(:, i)/(vars.ScaleSize * 2) + i * 2;
            YTickPos(i * 3 - 1) = i * 2 + 1;
            YTickPos(i * 3 - 2) = i * 2 + 0.5;
            YTickPos(i * 3) = i * 2 + 1.5;
            YTickLab(i * 3 - 1) = ChanStr(i);
            YTickLab(i * 3 - 2) = {sprintf('%.2f', -1 * vars.ScaleSize)};
            YTickLab(i * 3) = {sprintf('%.2f', vars.ScaleSize)};
        end
        
    end


    %make our graph
    if isgraphics(Graph.Graph)
        set(0, 'CurrentFigure', Graph.Graph)
        plot(time, SampleToPlot, 'k' )
        xLines = unique(round(xMin:xMax));
        hold on
        for xl = xLines
            xline(xl, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.25);
        end
        hold off
        xlim([xMin, xMax])
        ylim([1, vars.numChannelsToPlot * 2 + 2])
%         yticks(1 + (2:2:vars.numChannelsToPlot * 2));
%         if ~vars.UseKalman
%             yticklabels(vars.ChannelNames(vars.ChannelsToPlot))
%         end
        yticks(YTickPos);
        yticklabels(YTickLab);
        xlabel('Time (S)')
    end

    %set the next time we'll update our plot
    vars.EEGPlotPosition = vars.EEGPlotPosition + EEG.fs/vars.EEGPlotPerSecond;
end
end

