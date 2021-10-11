function [P] = ProbEyesClosedBayes(alpha, EODist, ECDist)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pEO = normpdf(alpha, EODist.mean, EODist.stdev);
pEC = normpdf(alpha, ECDist.mean, EODist.stdev);

P = pEC/(pEO+pEC);
end

