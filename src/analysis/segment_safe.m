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
% ignore_out_of_range : bool
%     If true, epochs out of range will be ignored (but wargning still issued).
% 
% Returns
% -------
% header_out : letswave6 header
% data_out : letswave6 data
% 

parser = inputParser;

addParameter(parser, 'x_start', []);
addParameter(parser, 'x_duration', []);
addParameter(parser, 'ignore_out_of_range', false);

parse(parser, varargin{:});

x_start = parser.Results.x_start;
x_duration = parser.Results.x_duration;
ignore_out_of_range = parser.Results.ignore_out_of_range;


if isempty(x_start)
    error('pls provide "x_start" as kwarg'); 
end
if isempty(x_duration)
    error('pls provide "x_duration" as kwarg'); 
end
%% 

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

% find matching event codes
codes_in_header = {events_in_header.code};
idx_matching_event = find(cellfun(@(x) ismember(x, segment_codes), ...
                          codes_in_header));

n_events = length(idx_matching_event);
n_x = round(x_duration / header.xstep);

data_out = nan(...
    n_events, size(data, 2), size(data, 3), size(data, 4), size(data, 5), n_x...
    ); 
                      
events_out = [];
c_event_out = 1;
c_epoch_out = 1;

out_of_range_events = [];

for i_event=1:n_events
      
    event_code = events_in_header(idx_matching_event(i_event)).code;
    event_epoch = events_in_header(idx_matching_event(i_event)).epoch;
    ep_latency = events_in_header(idx_matching_event(i_event)).latency + x_start;
    ep_onset_idx = round((ep_latency - header.xstart) / header.xstep); 
    
    if ep_onset_idx+n_x > header.datasize(end) || ep_onset_idx+1 < 1 
        if ignore_out_of_range
            warning(...
                'trying to segment from %.3f to %.3f but data is only from %.3f to %.3f ... ignoring this segment...', ...
                ep_latency, ...
                ep_latency + x_duration, ...
                header.xstart, ...
                header.datasize(end) * header.xstep - header.xstart ...
                );
            out_of_range_events = [out_of_range_events, i_event];
            continue
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
    end
    
    data_out(i_event, :, :, :, :, :) = ...
        data(event_epoch, :, :, :, :, ep_onset_idx+1 : ep_onset_idx+n_x);
    
    % check which other events fall into this new epoch
    other_events = events_in_header;
    other_events(idx_matching_event(i_event)) = [];
    
    other_events_idx = ...
        find([[other_events.epoch] == event_epoch & ...
             [other_events.latency] > ep_latency & ...
             [other_events.latency] < ep_latency + x_duration ...
             ]);
    
    other_events_idx_ignore = find(...
                      cellfun(@(x) ismember(x, segment_codes), ...
                      {other_events(other_events_idx).code})...
                      );
                  
    if ~isempty(other_events_idx_ignore)
        % ignore the events that are used for the current epoching        
          
        removed_codes = join(...
            {other_events(other_events_idx(other_events_idx_ignore)).code},...
            ', '...
            );
        
        warning('removing additional events "%s" from epoch %d', ...
                removed_codes{1}, i_event);

        other_events_idx(other_events_idx_ignore) = [];
                     
    end
    
    events_out(c_event_out).code = event_code;
    events_out(c_event_out).latency = 0;
    events_out(c_event_out).epoch = c_epoch_out;
    
    c_event_out = c_event_out + 1;
    
    % write additional events that fall into the segment
    for i_additional_event=1:length(other_events_idx)
        
        events_out(c_event_out).code = ...
            other_events(other_events_idx(i_additional_event)).code;
        
        events_out(c_event_out).latency = ...
            other_events(other_events_idx(i_additional_event)).latency - ...
            events_in_header(idx_matching_event(i_event)).latency;
        
        events_out(c_event_out).epoch = i_event;
        
        c_event_out = c_event_out + 1;
        
    end
    
    c_epoch_out = c_epoch_out + 1;
    
end

data_out(out_of_range_events, :, :, :, :, :) = [];

header_out = header;
header_out.datasize = size(data_out);
header_out.events = events_out;
header_out.xstart = x_start;



