function header = make_header(varargin)

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
for i=1:n_chans
    chanlocs(i).labels = sprintf('chan%d', i);
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