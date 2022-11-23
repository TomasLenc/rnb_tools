function plot_fft(freq, mX, varargin)

parser = inputParser; 

addParameter(parser, 'maxfreqlim', 5.1); 
addParameter(parser, 'frex_meter_rel', []); 
addParameter(parser, 'frex_meter_unrel', []); 
addParameter(parser, 'ax', []); 
addParameter(parser, 'linew', 1.7); 
addParameter(parser, 'fontsize', 12); 
addParameter(parser, 'col_meter_rel', [222 45 38]/255); 
addParameter(parser, 'col_meter_unrel', [49, 130, 189]/255); 
addParameter(parser, 'col_neutral', repmat(0.5, 1, 3)); 

parse(parser, varargin{:}); 

maxfreqlim = parser.Results.maxfreqlim; 
frex_meter_rel = parser.Results.frex_meter_rel; 
frex_meter_unrel = parser.Results.frex_meter_unrel; 
ax = parser.Results.ax; 
linew = parser.Results.linew; 
fontsize = parser.Results.fontsize; 
col_meter_rel = parser.Results.col_meter_rel; 
col_meter_unrel = parser.Results.col_meter_unrel; 
col_neutral = parser.Results.col_neutral; 


mX = ensure_row(mX); 
freq = ensure_row(freq); 
frex_meter_rel = ensure_row(frex_meter_rel); 
frex_meter_unrel = ensure_row(frex_meter_unrel); 

if isempty(ax)
    f = figure; 
    ax = axes(f); 
end
hold(ax, 'on');

stem(ax, freq, mX, 'color', col_neutral, 'marker', 'none', 'linew', linew); 

if ~isempty(frex_meter_rel)
    frex_meter_rel_idx = dsearchn(freq', frex_meter_rel')'; 
    
    stem(ax, freq(frex_meter_rel_idx), mX(frex_meter_rel_idx), ...
        'color', col_meter_rel,...
        'marker','none',...
        'linew',linew); 
end

if ~isempty(frex_meter_unrel)
    frex_meter_unrel_idx = dsearchn(freq', frex_meter_unrel')'; 
    
    stem(ax, freq(frex_meter_unrel_idx), mX(frex_meter_unrel_idx), ...
        'color', col_meter_unrel,...
        'marker','none',...
        'linew',linew); 
end

ax.Color = 'none'; 
ax.TickDir = 'out'; 
ax.XTick = []; 
ax.XLim = [0, maxfreqlim]; 
ax.YTick = [0, ax.YLim(2)]; 
ax.XAxis.Visible = 'off'; 
ax.FontSize = fontsize; 

