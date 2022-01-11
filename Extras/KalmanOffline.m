function [Recovered_eeg_Kalman] = KalmanOffline(EEG, K_Targets, q, R, varargin)
%   Performs Kalman filter offline to remove BCG artifact
% Inputs : EEG - data matrix, n samples x m + 2 channels (1 extra channel for
%               triggers, one for time. extra channels Can be filled with         
%               zeros if not present
%          K_Targets - a vector indicating which channel numbers to filter
%          q - a scalar hyperparameter, which must be tuned for optimal
%               performance. A reasonable starting point is 1.0e-4
%          R - a scalar hyperparameter, which must be tuned for optimal
%               performance. A reasonable starting point is 10000
%               Varargin - can be omitted, or set to a scalar to use a fixxed
%               number of reference channels. otherwise will autodetect the
%               correct number
% Outputs : Recovered_eeg_Kalman - clean data matrix, n samples x
%               length(K_targets) channel
                    

[n_all, nChans] = size(EEG); % detect dataset size
numEEG = length(K_Targets);
refChans = 1:nChans - 2;
auxChans = [25:32]'; % omit aux channels like ECG and EMG, E1, M1, from filter
refChans = setdiff(refChans, [K_Targets; auxChans]); % detect reference channels
if nargin > 4
    flag = 1;
else 
    flag = 0;
end
designMat_synth = [ones(n_all, 1), EEG]'; % build our design matrix from reference channels
DM = designMat_synth;
DM(2:end, :) = DM(2:end, :) - mean(DM(1+refChans, :), 1);
SynthEEG = EEG(:, K_Targets)';  % select target channel
% initialize algorithm variables
if flag
    j = size(varargin{1}, 2) + 1;
else
    j = length(refChans) + 1;
end
x_hat = zeros(j, numEEG);
h_hat = zeros(n_all, numEEG);
Q = eye(j) * q;
P = eye(j);
P = repmat(P, 1, 1, numEEG);
for chan = 1:numEEG % loop through target channels
    if  flag
        ChannelRefs = [1, 1 + varargin{1}(chan, :)];
        channelDM = DM(ChannelRefs, :);
    else
        channelDM = DM(refChans);
    end
    for t = 1:n_all % loop through time
        % perform Kalmfilter update on a sample - by - channel basis
        x_hat_givenPrevious = x_hat(:, chan);
        P_givenPrevious = P(:, :, chan) + Q;
        KalmanGain = P_givenPrevious * channelDM(:, t) ...
                    * (channelDM(:, t)' * P_givenPrevious * channelDM(:, t) + R) ^-1;
        x_hat(:, chan) = x_hat_givenPrevious + ...
                    KalmanGain * (SynthEEG(chan, t) - channelDM(:, t)' * x_hat_givenPrevious);
        P(:, :, chan) = (eye(j) - KalmanGain * channelDM(:, t)') * P_givenPrevious;
        h_hat(t, chan) = channelDM(:, t)' * x_hat(:, chan);
    end
end
% subtract modeled BCG from original EEG signal
Recovered_eeg_Kalman = (SynthEEG - h_hat')';
end


