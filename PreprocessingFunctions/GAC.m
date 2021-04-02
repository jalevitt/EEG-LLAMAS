function vars = GAC(EEG, vars)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if vars.UseGAC % Perform GAC
    if vars.trGap~= 0 || sum(EEG.Recording_orig(:, end - 1) == 253) > (vars.ntr + 1)
        if vars.trGap == 0
            trSamp = 1:(vars.currentPosition_orig - 1);
            trSamp = trSamp(EEG.Recording_orig(:, end - 1) == 253);
            vars.trGap = mode(diff(trSamp))
            EEG.trGap = vars.trGap;
        end
        %template = zeros(vars.trGap, ChansInChunk - 1,  vars.ntr);
        template = reshape(EEG.Recording_orig(vars.currentPosition_orig - vars.ntr * vars.trGap:vars.currentPosition_orig - 1, ...
            1:end - 2), vars.trGap, vars.ntr, vars.ChansInChunk - 1);
        template = permute(template, [1, 3, 2]);
        MeanTemplate = mean(template, 3);
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

