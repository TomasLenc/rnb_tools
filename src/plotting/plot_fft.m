function plot_fft(freq, mX, varargin)
% Plots magnitude spectra. 
%   
%   Parameters
%   ----------
%   freq :  array_like
%       Vector of frequencies in Hz. 
%   mX : array_like
%       Vector of FFT magniudes (one per frequency bin). 
%   maxfreqlim : float, default=5.1
%       Maximum frequency that will be plotted. 
%   frex_meter_rel : array_like, optional 
%       Meter/beat/whatever-related frequencies. 
%   frex_meter_unrel : array_like, optional 
%       Meter/beat/whatever-unrelated frequencies. 
%   ax : axis_object, optional 
%       Axis handle to plot in. 
%   linew : float, optional, default=1.7
%       Line width. 
%   fontsize : float, optional, default=12
%       Font sizwe. 
%   col_meter_rel : [float, float, float], optional, default=[222 45 38]/255
%       RGB triplet to use as a color for meter-related frequwncies. 
%   col_meter_unrel : [float, float, float], optional, default=[49, 130, 189]/255
%       RGB triplet to use as a color for meter-unrelated frequwncies. 
%   col_neutral : [float, float, float], optional, default=[.5 .5 .5]
%       RGB triplet to use as a baseline color for all frequencies. 
%   
%   Returns
%   -------
%   -

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

