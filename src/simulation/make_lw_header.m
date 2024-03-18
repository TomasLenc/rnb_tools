function header = make_lw_header(varargin)

parser = inputParser;

addParameter(parser, 'data', []);
addParameter(parser, 'x_start', 0);
addParameter(parser, 'fs', nan);
addParameter(parser, 'n_chans', nan);
addParameter(parser, 'events', []);

parse(parser, varargin{:});

data = parser.Results.data;
x_start = parser.Results.x_start;
fs = parser.Results.fs;
n_chans = parser.Results.n_chans;
events = parser.Results.events;

chanlocs = [];
if ~isempty(data) && isnan(n_chans)
    n_chans = size(data, 2); 
end
for i=1:n_chans
    chanlocs(i).labels = sprintf('chan%d', i);
    chanlocs(i).topo_enabled = false; 
    chanlocs(i).SEEG_enabled = false; 
end

if isempty(events)
    for i_trial=1:size(data, 1)
        events(i_trial).code = 'trial';
        events(i_trial).latency = 0;
        events(i_trial).epoch = i_trial;
    end
end

header = [];

header.filetype = 'time_amplitude';
header.name = 'test';
header.tags = '';
header.history = [];
header.datasize = size(data);
header.xstart = x_start;
header.ystart = 0;
header.zstart = 0;
header.xstep = 1/fs;
header.ystep = 1;
header.zstep = 1;
header.chanlocs = chanlocs;
header.events = events;
