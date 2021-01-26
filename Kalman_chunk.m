function [Kalman, P, x_hat] = Kalman_chunk(chunk, P, Q, R, r, bcg_idx, x_hat, KalmanTargets)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
chunk = chunk - mean(chunk, 1); % re-reference to avereage
numReg = length(r);
n = size(chunk, 2);
BCG_all = chunk(bcg_idx, :);
BCG = BCG_all(r, :);
DM = ones(numReg + 1, n);
DM(2:end, :) = BCG;
numKalman = length(KalmanTargets);
Signal = chunk(KalmanTargets, :);
h_hat = zeros(n, numKalman);
for chan = 1:numKalman
    for t = 1:n
        x_hat_givenPrevious = x_hat(:, chan);
        P_givenPrevious = P(:, :, chan) + Q;
        KalmanGain = P_givenPrevious * DM(:, t) ...
            * (DM(:, t)' * P_givenPrevious * DM(:, t) + R)^-1;
        x_hat(:, chan) = x_hat_givenPrevious + ...
            KalmanGain * (Signal(chan, t) - DM(:, t)' * x_hat_givenPrevious);
        P(:, :, chan) = (eye(numReg + 1) - KalmanGain * DM(:, t)') * P_givenPrevious;
        h_hat(t, chan) = DM(:, t)' * x_hat(:, chan);
    end
    
end
Kalman = Signal - h_hat';
end

