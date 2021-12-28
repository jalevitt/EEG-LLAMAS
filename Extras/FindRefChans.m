function RefMatrix = FindRefChans(targets, refs, numRefs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load('ChannelLocs.mat', 'locs')

targetLocs = [[locs(targets).X]; [locs(targets).Y]; [locs(targets).Z]]';
RefLocs = [[locs(refs).X]; [locs(refs).Y]; [locs(refs).Z]]';

RefMatrix = zeros(length(targets), numRefs);

for c = 1:length(targets)
   currentLoc = targetLocs(c, :);
   DistanceToTargets = vecnorm(RefLocs - currentLoc, 2, 2);
   [~, FarthestChans] = maxk(DistanceToTargets, numRefs);
   RefMatrix(c, :) = refs(FarthestChans);
end

end

