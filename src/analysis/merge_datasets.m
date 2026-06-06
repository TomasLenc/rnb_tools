function [header, data] = merge_datasets(headers, datas)
% MERGE_DATASETS  Merge multiple Letswave6 datasets into a single dataset.
%
% This function concatenates multiple Letswave6 datasets along the 6th
% dimension and merges their event structures into a unified header.
% Event latencies from subsequent datasets are shifted to reflect their
% position in the concatenated signal.
%
% Parameters
% ----------
% headers : cell array of struct
%     Cell array of Letswave6 header structures. Each header must include:
%     - datasize (numeric array)
%     - xstep (double)
%     - chanlocs (struct array with field 'labels')
%     - events (struct array, optional)
%
% datas : cell array of numeric arrays
%     Cell array of data arrays corresponding to each header. All arrays
%     must be compatible for concatenation along the 6th dimension.
%
% Returns
% -------
% header : struct
%     Merged header. Event latencies are adjusted based on cumulative
%     signal duration. Channel information is preserved from the first dataset.
%
% data : numeric array
%     Concatenated data array along the 6th dimension.
%
% Notes
% -----
% - headers and datas must have the same length.
% - Channel labels must match exactly (including order) across datasets.
% - Event structs are copied fully; all additional fields are preserved.
% - Event latencies are shifted by:
%       previous_datasize(end) * xstep
% - All merged events are assigned to epoch = 1.
% - Assumes all datasets differ only along the 6th dimension.

% basic checks
if numel(headers) ~= numel(datas)
    error('merge_datasets:InputSizeMismatch', ...
        'headers and datas must have the same number of elements.');
end

header = headers{1};
data   = datas{1};

% ensure events field exists
if ~isfield(header, 'events') || isempty(header.events)
    header.events = [];
end

% reference channel labels
labels = {header.chanlocs.labels};

for i_dataset = 2:numel(datas)

    current_header = headers{i_dataset};

    % Check channel consistency
    if numel(current_header.chanlocs) ~= numel(labels)
        error('merge_datasets:ChannelCountMismatch', ...
            'Channel count mismatch between datasets.');
    end

    current_labels = {current_header.chanlocs.labels};

    if ~isequal(labels, current_labels)
        error('merge_datasets:ChannelLabelMismatch', ...
            'Channel labels mismatch between datasets.');
    end

    % Merge events
    if isfield(current_header, 'events') && ~isempty(current_header.events)

        for i_event = 1:numel(current_header.events)

            new_event = current_header.events(i_event);

            new_event.latency = ...
                (header.datasize(end) * header.xstep) + ...
                current_header.events(i_event).latency;

            new_event.epoch = 1;

            header.events(end + 1) = new_event; %#ok<AGROW>
        end
    end

    % Concatenate data
    data = cat(6, data, datas{i_dataset});

    % Update datasize
    header.datasize(end) = header.datasize(end) + ...
                           current_header.datasize(end);

end