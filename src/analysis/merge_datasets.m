function [header, data] = merge_datasets(datasets) 
% MERGE_DATASETS Merge multiple Letswave6 datasets into a single dataset.
%
% This function concatenates the data arrays from multiple Letswave6
% datasets along the 6th dimension and combines their event structures
% into a unified header. Event latencies from subsequent datasets are
% shifted to reflect their new positions in the merged data.
%
% Parameters
% ----------
% datasets : struct array
%     Array of Letswave6 dataset structures. Each element must contain:
%     - header : struct
%         Metadata describing the dataset, including:
%         - datasize (numeric array): Size of the data dimensions.
%         - xstep (double): Temporal resolution (sampling interval).
%         - events (struct array): Event information with fields:
%             - latency (double): Event time position.
%             - code (char or numeric): Event identifier.
%     - data : numeric array
%         Multidimensional data array. All datasets must be compatible
%         for concatenation along the 6th dimension.
%
% Returns
% -------
% header : struct
%     Merged header structure. Event latencies are adjusted such that
%     events from subsequent datasets are shifted by the cumulative
%     duration of preceding datasets. The datasize field is updated
%     accordingly.
%
% data : numeric array
%     Merged data array created by concatenating input datasets along
%     the 6th dimension.
%
% Notes
% -----
% - The first dataset is used as the base dataset.
% - Event latencies are shifted by:
%       previous_datasize(end) * xstep
% - All merged events are assigned to epoch = 1.
% - Assumes datasets only differ along the 6th dimension.
%

header = datasets(1).header; 
data = datasets(1).data; 

labels = {header.chanlocs.labels};

for i_datset=2:length(datasets)

    current_header = datasets(i_datset).header; 
    
    % ensure that channel labels match
    if ~isequal(labels, {current_header.chanlocs.labels})
        error('merge_datasets:ChannelLabelMismatch', ...
            'Channel labels mismatch between datasets.');
    end

    % Merge events
    for i_event = 1:length(current_header.events)

        % copy full event struct
        new_event = current_header.events(i_event);

        % adjust latency
        new_event.latency = ...
            (header.datasize(end) * header.xstep) + ...
            current_header.events(i_event).latency;

        % enforce merged dataset properties
        new_event.epoch = 1;

        % append event
        header.events(end + 1) = new_event;

    end

    data = cat(6, data, datasets(i_datset).data); 
    
    header.datasize(end) = header.datasize(end) + ...
                                datasets(i_datset).header.datasize(end); 

end
