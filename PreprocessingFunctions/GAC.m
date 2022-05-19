function vars = GAC(EEG, vars)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if vars.UseGAC % Perform GAC
    
    if  sum(EEG.Recording_orig(:, end - 1) == vars.trMarker) > (vars.ntr + 1) % || vars.trGap~= 0 
        if vars.trGap == 0 % if we haven't yet detected the samples per TR, do so now, store in vars.trGap
            trSamp = 1:(vars.currentPosition_orig - 1);
            trSamp = trSamp(EEG.Recording_orig(:, end - 1) == vars.trMarker);
            vars.trGap = mode(diff(trSamp))
            EEG.trGap = vars.trGap;
        end
        
        % get the data from the last n TRs, and reshape so its vars.trGap
        % samples X M channels X vars.ntr epochs
        template = reshape(EEG.Recording_orig(vars.currentPosition_orig - vars.ntr * vars.trGap:vars.currentPosition_orig - 1, ...
            1:end - 2), vars.trGap, vars.ntr, vars.ChansInChunk - 1); %
        template = permute(template, [1, 3, 2]);
        
        % get template average across n TRs
        MeanTemplate = mean(template, 3);
        
        % basically just subtract the template from our data - this
        % implemnation works in real time, your implementation will differ
        % offline
        if vars.SamplesInChunk <= vars.trGap
            vars.OrigChunk(1:end - 1, :) = vars.OrigChunk(1:end - 1, :) - ...
                MeanTemplate(end - vars.SamplesInChunk + 1:end, :)';
        else
            numReps = floor(vars.SamplesInChunk / vars.trGap) + 1;
            MeanTemplate = repmat(MeanTemplate, numReps, 1);
            vars.OrigChunk(1:end - 1, :) = vars.OrigChunk(1:end - 1, :) - ...
                MeanTemplate(end - vars.SamplesInChunk + 1:end, :)';

        end
    end
end
end

