function [vars, Graph, EEG] = Kalman_chunk(EEG,vars, Graph)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if vars.UseKalman
    chunk = vars.chunk - mean(vars.chunk, 1); % re-reference to avereage
    numReg = length(vars.r);
    n = size(chunk, 2);
    BCG = chunk(vars.r, :);
    %BCG = BCG_all(vars.r, :);
    DM = ones(numReg + 1, n);
    DM(2:end, :) = BCG;
    numKalman = length(vars.KalmanTargets);
    Signal = chunk(vars.KalmanTargets, :);
    h_hat = zeros(n, numKalman);
    for chan = 1:numKalman
        for t = 1:n
            x_hat_givenPrevious = vars.x_hat(:, chan);
            P_givenPrevious = vars.P(:, :, chan) + vars.Q;
            KalmanGain = P_givenPrevious * DM(:, t) ...
                * (DM(:, t)' * P_givenPrevious * DM(:, t) + vars.R)^-1;
            vars.x_hat(:, chan) = x_hat_givenPrevious + ...
                KalmanGain * (Signal(chan, t) - DM(:, t)' * x_hat_givenPrevious);
            vars.P(:, :, chan) = (eye(numReg + 1) - KalmanGain * DM(:, t)') * P_givenPrevious;
            h_hat(t, chan) = DM(:, t)' * vars.x_hat(:, chan);
        end

    end
    Kalman = Signal - h_hat';
    EEG.Kalman_Signal(vars.currentPosition - vars.SamplesInChunk:vars.currentPosition - 1, :) ...
                            = Kalman';
end
    
end
