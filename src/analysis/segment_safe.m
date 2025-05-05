function [header_out, data_out] = segment_safe(header, data, segment_codes, varargin)
% This is a replacement for the unsafe RLW_segmentation function from
% letswave6. When asked to segment out of data range, RLW_segmentation doesn't
% even issue a warning, and instead returns an epoch full of zeros. Here, we
% explicitly raise an exception in such situation, unless a special flag is
% passed, in which case we still issue a warning. 
%
% Parameters
% ----------
% header : letswave6 header
% data : letswave6 data
% segment_codes : cell of strings
%     Event codes that will be used for segmentation. 
% out_of_range : str, {'error', 'ignore', 'zero_pad'}, optional, default='error' 
%     If the requested epoch goes out of the range of the continous data,
%     we need to decide what to do (original lw6 simply returned zeroes for
%     that epoch without any warning). Here, we can either let the function
%     raise an error, ignore the out-of-range epoch, or pad the
%     out-of-range portion of that epoch with zeroes. 
% keep_other_events : bool, optional, default=true
%     If true, any events with a code that's not in `segment_codes` will be
%     kept in the epoched data. This can be useful to, e.g. keep an event
%     about behavioral response in each trial etc. 
% 
% Returns
% -------
% header_out : letswave6 header
% data_out : letswave6 data
% 

parser = inputParser;

addParameter(parser, 'x_start', []);
addParameter(parser, 'x_duration', []);
addParameter(parser, 'out_of_range', 'error'); % 'error', 'ignore', 'zero_pad'
addParameter(parser, 'keep_other_events', true); 

parse(parser, varargin{:});

x_start = parser.Results.x_start;
x_duration = parser.Results.x_duration;
out_of_range_treatment = parser.Results.out_of_range;
keep_other_events = parser.Results.keep_other_events;


if isempty(x_start)
    error('pls provide "x_start" as kwarg'); 
end
if isempty(x_duration)
    error('pls provide "x_duration" as kwarg'); 
end
%% 

fprintf('\nchecking for event duplicates... ');

% only look at unique events
events_in_header = header.events(1);
for i_event=1:length(header.events)
    if ~(...
        ismember(header.events(i_event).code, {events_in_header.code}) && ...
        ismember(header.events(i_event).latency, [events_in_header.latency]) && ...
        ismember(header.events(i_event).epoch, [events_in_header.epoch])...
        )
        events_in_header(end+1) = header.events(i_event); 
    end
end
fprintf('%d duplicated events found...\n', numel(events_in_header) - numel(header.events));

% find event codes matching the ones requested for segmentation by the user
codes_in_header = {events_in_header.code};

idx_matching_event = find(cellfun(@(x) ismember(x, segment_codes), ...
                          codes_in_header));

n_ep = length(idx_matching_event);

% number of samples for the segment
n_x = round(x_duration / header.xstep);

% allocate data output with correct shape
data_out = nan(...
    n_ep, size(data, 2), size(data, 3), size(data, 4), size(data, 5), n_x...
    ); 
       
% allocate output event structure  
event_fieldnames = fieldnames(events_in_header); 
event_fieldnames = [event_fieldnames; {'latency_orig'}]; 

% inti counter for output event sturct
c_event_out = 1;

% also we need another counter that will keep track of which epoch in the
% output array we are currently at. Note that this can be different from
% the event number as some events may be skipped due to being out-of-range.
% Likewise, we may want to keep other events (e.g. response events) from the
% continuous data - so one epoch can end up having several triggers inside
% it. 
c_epoch_out = 1;

% Here we will save the indices of events that ended up being out of range
% of the original data array. Later, we can use this to remove those epochs
% from the output data array (as there will be nothing assigned, just nans
% that we allocated)
out_of_range_events = [];

fprintf('preparing epoched data...\n');

% go over each event that contains a code matching what the user wants to
% use for segmentation  
for i_event=1:n_ep
      
    % take that matching event out of the header
    event = events_in_header(idx_matching_event(i_event)); 

    % check the epoch where this event is in the input data array
    event_epoch = event.epoch;
    
    % the latency where segmentation will start 
    ep_latency = event.latency + x_start;
    
    % get the sample index for the begining of the segment 
    ep_onset_idx = round((ep_latency - header.xstart) / header.xstep); 
    
    % check if any part of the segment is out of data range and issue warning/error
    if ep_onset_idx+n_x > header.datasize(end) || ep_onset_idx+1 < 1 
        
        if strcmp(out_of_range_treatment, 'ignore')
            
            warning(...
                'trying to segment from %.3f to %.3f but data is only from %.3f to %.3f ... ignoring this segment...', ...
                ep_latency, ...
                ep_latency + x_duration, ...
                header.xstart, ...
                header.datasize(end) * header.xstep - header.xstart ...
                );
            out_of_range_events = [out_of_range_events, i_event];
            continue
            
        elseif strcmp(out_of_range_treatment, 'zero_pad')
            
            warning(...
                'trying to segment from %.3f to %.3f but data is only from %.3f to %.3f ... \nI will zero pad the missing data ...', ...
                ep_latency, ...
                ep_latency + x_duration, ...
                header.xstart, ...
                header.datasize(end) * header.xstep - header.xstart ...
                );
            
            % prepare zero padding 
            if ep_onset_idx+1 < 1 
                n_pad = - ep_onset_idx; 
                data_start_idx = 1; 
                pad_shape = header.datasize; 
                pad_shape(6) = n_pad; 
                pad_before = zeros(pad_shape);                
            else
                pad_before = []; 
                data_start_idx = ep_onset_idx+1; 
            end
            
            if ep_onset_idx+n_x > header.datasize(end)
                n_pad = (ep_onset_idx+n_x) - header.datasize(end); 
                data_end_idx = header.datasize(end); 
                pad_shape = header.datasize; 
                pad_shape(6) = n_pad; 
                pad_after = zeros(pad_shape); 
            else
                pad_after = []; 
                data_end_idx = ep_onset_idx+n_x; 
            end
            
            % extract the data setment, zero pad, and assign to the output array
            data_out(i_event, :, :, :, :, :) = ...
                cat(6, pad_before, ...
                       data(event_epoch, :, :, :, :, data_start_idx : data_end_idx), ...
                       pad_after);
            
        else
            
            error(...
                'segment_safe:SegmentOutOfDataRange', ...
                'trying to segment from %.3f to %.3f but data is only from %.3f to %.3f', ...
                ep_latency, ...
                ep_latency + x_duration, ...
                header.xstart, ...
                header.datasize(end) * header.xstep - header.xstart ...
                );
        end
        
    else
    
        % assign the data setment to the output array
        data_out(i_event, :, :, :, :, :) = ...
            data(event_epoch, :, :, :, :, ep_onset_idx+1 : ep_onset_idx+n_x);
        
    end
    
    % --------------------------------------------------------------------
    % write event to output structure 
    
    for i_name=1:length(event_fieldnames)
        if isfield(event, event_fieldnames{i_name})
            events_out(c_event_out).(event_fieldnames{i_name}) = event.(event_fieldnames{i_name}); 
        end
    end
    
    % fix latency to 0
    events_out(c_event_out).latency = 0; 
    
    % fix epoch number to the correct epoch in the output data array
    events_out(c_event_out).epoch = c_epoch_out;
    
    % include original latency 
    events_out(c_event_out).latency_orig = event.latency; 
    
    % update the event counter
    c_event_out = c_event_out + 1;
    
    % --------------------------------------------------------------------
    % other events 
    
    if keep_other_events
        
        % check if there are any other events falling into this new epoch
        other_events = events_in_header;
        other_events(idx_matching_event(i_event)) = [];

        other_events_idx = ...
            find([[other_events.epoch] == event_epoch & ...
                 [other_events.latency] > ep_latency & ...
                 [other_events.latency] < ep_latency + x_duration ...
                 ]);

        % It may happen that the next/prev event that the user is asking to use
        % for segmentation falls into the current epoch (this can be for
        % example when one is segmenting trials with a buffer time before/after
        % the trial, and there is prev/next trial event that falls into that
        % time range). Let's get rid of those events. 
        other_events_idx_ignore = find(...
                          cellfun(@(x) ismember(x, segment_codes), ...
                          {other_events(other_events_idx).code})...
                          );

        % ignore the events that are used for the current epoching (and issue
        % warning)
        if ~isempty(other_events_idx_ignore)

            removed_codes = join(...
                {other_events(other_events_idx(other_events_idx_ignore)).code},...
                ', '...
                );

            warning('removing additional events "%s" from epoch %d', ...
                    removed_codes{1}, i_event);

            other_events_idx(other_events_idx_ignore) = [];

        end

        % write additional events that fall into the segment
        for i_add_event=1:length(other_events_idx)

            event = other_events(other_events_idx(i_add_event)); 
            
            for i_name=1:length(event_fieldnames)
                if isfield(event, event_fieldnames{i_name})
                    events_out(c_event_out).(event_fieldnames{i_name}) = event.(event_fieldnames{i_name}); 
                end
            end
           
            % fix latency so it is wrt the event used for segmentation 
            events_out(c_event_out).latency = ...
                    events_out(c_event_out).latency - ...
                    events_in_header(idx_matching_event(i_event)).latency;

                
            % include original latency 
            events_out(c_event_out).latency_orig = event.latency; 
                
            % set the correct epoch in the output array
            events_out(c_event_out).epoch = i_event;

            % update the counter for output event structure 
            c_event_out = c_event_out + 1;

        end
    
    end
    
    % --------------------------------------------------------------------
    % update epoch counter 
    c_epoch_out = c_epoch_out + 1;
        
end

% remove slices that correcpond to events which ended up giving out of
% range epochs (there are just nans there in the output array anyway)
data_out(out_of_range_events, :, :, :, :, :) = [];

% put the output header together
header_out = header;
header_out.datasize = size(data_out);
header_out.events = events_out;
header_out.xstart = x_start;

fprintf('done...\n');


