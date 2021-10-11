function occchan=occ_chan_per_subj(subj)
% selected manually as a channel that has high occipital alpha at baseline,
% and looks high quality in cleaned MR data
% run allchans_eyesOC

switch subj
    case 'osceeg2b'
        occchan='c126'; % not great
    case 'osceeg3b'
        occchan='c136';
    case 'osceeg5b'
        occchan='c139';
    case 'osceeg6b'
        occchan='c136';
    case 'osceeg7b'
        occchan='c116';
    case 'osceeg8b'
        occchan='c126';
    case 'osceeg9b'
        occchan='c124';
    case 'osceeg10b'
        occchan='c124';
    case 'osceeg15b'
        occchan='c136';
    case 'osceeg17b'
        occchan='c118';
    case 'osceeg19b'
        occchan='c138';
    case 'osceeg20b'
        occchan='c125'; % look into quality further
    case 'osceeg21b'
        occchan='c137';
    case 'osceeg22b'
        occchan='c116';
    case 'inkrefcap1'
        occchan='c116';
    case 'inkrefcap2'
        occchan='c126';
    case 'inkrefcap3'
        occchan='c126'; % pretty bad
    case 'inkrefcap4'
        occchan='c150'; % bad
    case 'inkrefcap5'
        occchan='c117';
    case 'inkrefcap6'
        occchan='c116';
    case 'inkrefcap7b'
        occchan='c118';
    case 'inkrefcap8'
        occchan='c118';
    case 'inkrefcap9a'
        occchan='c126';
end