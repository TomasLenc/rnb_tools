function [header, data] = interpolate_bad_chans(header, data, bad_chans, n_chans_interp)

good_chans_mask = cellfun(@(x) ~any(strcmpi(x, bad_chans)), ...
                        {header.chanlocs.labels}); 

good_chans = header.chanlocs(good_chans_mask); 

for i_bad=1:length(bad_chans)

    bad_chan_mask = strcmpi({header.chanlocs.labels}, bad_chans{i_bad});

    bad_chan = header.chanlocs(bad_chan_mask); 

    dist = nan(1, length(good_chans));        

    for i=1:length(good_chans)
        dist(i) = sqrt((good_chans(i).X - bad_chan.X)^2 ...
                       + (good_chans(i).Y - bad_chan.Y)^2 ...
                       + (good_chans(i).Z - bad_chan.Z)^2);
    end

    [~, idx_closest] = sort(dist);
    closest_channels = good_chans(idx_closest(1 : n_chans_interp));

    % if any of the closest channels are within the bad channel
    % list, issue a warning
    closest_match_bad_mask = cellfun(@(x) any(strcmpi(x, bad_chans)),...
                                     {closest_channels.labels}); 

    if any(closest_match_bad_mask)
        error('Interpolating ch %s: closest ch %s overlaps with another bad channel...', ...
               bad_chan.labels, closest_channels(closest_match_bad_mask).labels); 
    end

    tmp = join({closest_channels.labels}, ' '); 
    fprintf('Interpolating ch %s with avg of %s\n', bad_chan.labels,  tmp{1}); 

    [header, data] = RLW_interpolate_channel(header, data, ...
                            bad_chan.labels, {closest_channels.labels}); 

end