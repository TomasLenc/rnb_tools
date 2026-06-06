function [out_header, out_data] = crop_safe(header, data, varargin)
% CROP_SAFE  Crop a Letswave6 dataset along the time (x) dimension.
%
% This function extracts a time window from a Letswave6 dataset and
% optionally resets the time axis so that the cropped data starts at t = 0.
%
% Parameters
% ----------
% header : struct
%     Letswave6 header structure (must include datasize, xstart, xstep, events).
%
% data : numeric array
%     Multidimensional data array.
%
% varargin : key/value pairs
%     'x_start'       (double)  Start time in seconds (default = 0)
%     'x_end'         (double)  End time in seconds (default = end of data)
%     'reset_xstart'  (logical) If true, new xstart = 0 (default = true)
%
% Returns
% -------
% out_header : struct
%     Updated header with adjusted datasize, xstart, and events.
%
% out_data : numeric array
%     Cropped data array.
%
% Notes
% -----
% - Events outside the selected window are removed.
% - If reset_xstart = true:
%       new_latency = old_latency - x_start
% - If reset_xstart = false:
%       event latencies remain unchanged.

% Input parsing and setup
parser = inputParser;
parser.FunctionName = 'crop_safe';

default_x_start = 0;
default_x_end = header.xstart + (header.datasize(6) - 1) * header.xstep;
default_reset = true;

addParameter(parser, 'x_start', default_x_start, @isnumeric);
addParameter(parser, 'x_end', default_x_end, @isnumeric);
addParameter(parser, 'reset_xstart', default_reset, @islogical);

parse(parser, varargin{:});

x_start = parser.Results.x_start;
x_end   = parser.Results.x_end;
reset   = parser.Results.reset_xstart;

%%

% convert time → indices
dx1 = round((x_start - header.xstart) / header.xstep) + 1;
dx2 = round((x_end   - header.xstart) / header.xstep) + 1;

% clamp to valid time range
dx1 = max(1, dx1);
dx2 = min(header.datasize(6), dx2);

% crop data
out_data = data(:, :, :, :, :, dx1:dx2);

% update header
out_header = header;
out_header.datasize(6) = dx2 - dx1 + 1;

if reset
    out_header.xstart = 0;
else
    out_header.xstart = x_start;
end

% adjust events
out_events = [];
k = 1;

if isfield(header, 'events') && ~isempty(header.events)
    for i = 1:numel(header.events)
        old_latency = header.events(i).latency;

        if old_latency >= x_start && old_latency <= x_end
            out_events{k} = header.events(i); %#ok<AGROW>

            if reset
                out_events{k}.latency = old_latency - x_start;
            else
                out_events{k}.latency = old_latency;
            end

            k = k + 1;
        end
    end
end

if isempty(out_events)
    warning('no events inside the cropped time range')
    out_header.events = []; 
else
    out_header.events = cat(2, out_events{:});
end

% clean header
if isfield(out_header, 'epochdata')
    out_header = rmfield(out_header, 'epochdata');
end
