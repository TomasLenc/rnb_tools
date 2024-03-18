function idx = get_chan_idx(header, chan_labels)

idx = cellfun(@(x)find(strcmpi({header.chanlocs.labels},x)),...
                           chan_labels, 'uni', 1);
