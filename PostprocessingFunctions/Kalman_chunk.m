function [vars, Graph, EEG] = Kalman_chunk(EEG,vars, Graph)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    chunk = vars.chunk;
    chunk(vars.KalmanTargets, :) = chunk(vars.KalmanTargets, :) - mean(chunk(vars.KalmanTargets, :), 1); % re-reference EEG channels of chunk to average
    numReg = vars.numReg;
    n = size(chunk, 2);
    
    % build Design Matrix
    DM = ones(size(chunk, 1) + 1, n);
    DM(2:end, :) = chunk - mean(chunk(vars.r, :), 1);
    numKalman = length(vars.KalmanTargets);
    Signal = chunk(vars.KalmanTargets, :);
    h_hat = zeros(n, numKalman);
    
    % loop through channels and samples
    for chan = 1:numKalman
        if  strcmp(vars.RefMode, 'Channel Location')
            ChannelRefs = [1, 1 + vars.RefChans(chan, :)];
            channelDM = DM(ChannelRefs, :);
            
        else
            channelDM = DM([1; vars.r + 1], :);
        end
        for t = 1:n
            x_hat_givenPrevious = vars.x_hat(:, chan);
            P_givenPrevious = vars.P(:, :, chan) + vars.Q;
            KalmanGain = P_givenPrevious * channelDM(:, t) ...
                * (channelDM(:, t)' * P_givenPrevious * channelDM(:, t) + vars.R)^-1;
            vars.x_hat(:, chan) = x_hat_givenPrevious + ...
                KalmanGain * (Signal(chan, t) - channelDM(:, t)' * x_hat_givenPrevious);
            vars.P(:, :, chan) = (eye(numReg + 1) - KalmanGain * channelDM(:, t)') * P_givenPrevious;
            h_hat(t, chan) = channelDM(:, t)' * vars.x_hat(:, chan);
        end

    end
    Kalman = Signal - h_hat';
    EEG.Kalman_Signal(vars.currentPosition - vars.SamplesInChunk:vars.currentPosition - 1, :) ...
                            = Kalman';
% end
    
end
