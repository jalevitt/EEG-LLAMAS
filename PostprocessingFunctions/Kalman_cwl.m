function [vars, Graph, EEG] = Kalman_cwl(EEG,vars, Graph)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    if ~isfield(vars, 'useCWL')
        vars.useCWL = 1;
        vars.r = 33:36;
        vars.numReg = length(vars.r);
        vars.x_hat = zeros(vars.numReg + 1, length(vars.KalmanTargets));
        vars.Q = eye(vars.numReg + 1) * vars.Q(1, 1);
        vars.P = eye(vars.numReg + 1);
        vars.P = repmat(vars.P, 1, 1, length(vars.KalmanTargets));
        EEG.KalmanRegressors = vars.r;
        vars.lags = -5:0; % currently this will cause problems if max(lags) > 0
    end
    
    chunkSize = size(vars.chunk, 2) + abs(min(vars.lags));
    chunk = vars.chunk(:, (vars.currentPostion - chunkSize):(vars.currentPosition - 1));
    % chunk(vars.KalmanTargets, :) = chunk(vars.KalmanTargets, :) - mean(chunk(vars.KalmanTargets, :), 1); % re-reference EEG channels of chunk to average
    vars.numReg = 4;
    numReg = vars.numReg;
    n = size(chunk, 2);
    vars.cwl_chans = 33:36;
    
    % build Design Matrix
    DM = chunk(:, vars.cwl_chans);
    numKalman = length(vars.KalmanTargets);
    Signal = chunk(vars.KalmanTargets, :);
    h_hat = zeros(n, numKalman);
    
    % loop through channels and samples
    for chan = 1:numKalman
        channelDM = DM;
        for t = (1 + abs(min(vars.lags))):(n - abs(max(vars.lags)))
            temp = channelDM(:, t + vars.lags);
            temp = [1; temp(:)];
            x_hat_givenPrevious = vars.x_hat(:, chan);
            P_givenPrevious = vars.P(:, :, chan) + vars.Q;
            KalmanGain = P_givenPrevious * temp ...
                * (temp' * P_givenPrevious * temp + vars.R)^-1;
            vars.x_hat(:, chan) = x_hat_givenPrevious + ...
                KalmanGain * (Signal(chan, t) - temp' * x_hat_givenPrevious);
            vars.P(:, :, chan) = (eye(numReg + 1) - KalmanGain * temp') * P_givenPrevious;
            h_hat(t, chan) = temp' * vars.x_hat(:, chan);
        end

    end
    Kalman = Signal - h_hat';
    Kalman = Kalman(:, (1 + abs(min(vars.lags))):(n - abs(max(vars.lags))));
    Kalman = Kalman - mean(Kalman, 1);
    EEG.Kalman_Signal(vars.currentPosition - vars.SamplesInChunk:vars.currentPosition - 1, :) ...
                            = Kalman';
% end
    
end
