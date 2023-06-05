function h = plot_erp(x, varargin)
% This function plots time-domain ERP. Note that atm, this only plots 1-D responses, 
% not 2-D time-frequency heatmaps. 
% 
% If multiple epochs are passed, this will plot the average +- bootstrapped
% CIs as shaded regions. Otherwise it will just plot a line. 
% 
% Parameters
% ----------
% x : array 
%   (epoch or subject, time) array with the time-domain data
% fs : float
%   sampling rate 
% y_max : array | float | str
%   Either two element array with minimum and maximum y-axis limit. If only one
%   float is passed, it has to be positive, and then the axis limits will be
%   set from - to + that value. If 'any(isinf(yLim))', limits will be set based on the
%   data. 
% col : array 
%   RGB triplet for line color 
% varargin
%   'no_ci' : arg, str
%       don't plot bootstrapped confidence intervals when passing multiple
%       epochs or subjects in the input 


parser = inputParser; 

addParameter(parser, 'ax', []); 
addParameter(parser, 'fs', []); 
addParameter(parser, 't', []); 
addParameter(parser, 'y_max', Inf); 
addParameter(parser, 'sound_fs', []); 
addParameter(parser, 'sound_s', []); 
addParameter(parser, 'win_start_times', []); 
addParameter(parser, 'win_end_times', []); 
addParameter(parser, 'linew', 2); 
addParameter(parser, 'col', [0, 0, 0]); 
addParameter(parser, 'col_sound', [0.7, 0.7, 0.7]); 
addParameter(parser, 'col_win', [240, 174, 93]/255); 
addParameter(parser, 'ci', false); 
addParameter(parser, 'ci_method', 'analytical'); % analytical, boot
addParameter(parser, 'fontsize', 12); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
fs = parser.Results.fs; 
t = parser.Results.t; 
y_max = parser.Results.y_max; 
sound_s = parser.Results.sound_s; 
sound_fs = parser.Results.sound_fs; 
win_start_times = parser.Results.win_start_times; 
win_end_times = parser.Results.win_end_times; 
col = parser.Results.col; 
col_sound = parser.Results.col_sound; 
col_win = parser.Results.col_win; 
linew = parser.Results.linew; 
plot_ci = parser.Results.ci; 
ci_method = parser.Results.ci_method; 
fontsize = parser.Results.fontsize; 

if isempty(fs) && isempty(t)
    error('You must specify either "fs" or "t" parameter.'); 
end


%%

if isempty(ax)
    f = figure; 
    ax = axes(f); 
end

hold(ax,'on');

% if x is 1D, make sure it's a row vector
x = ensure_row(x);

if isempty(t)
    t = [0 : size(x, 2) - 1] / fs;
end

% average epochs (if any)
x_avg = mean(x, 1); 

% plot sound
if ~isempty(sound_fs) && ~isempty(sound_s)
    t_s = [0 : length(sound_s)-1] / sound_fs; 
    if isnumeric(y_max) && ~isinf(y_max)
     sound_s = sound_s / max(sound_s) * y_max;    
    else
        x_range = max(x_avg) - min(x_avg);
        sound_s = sound_s - min(sound_s); 
        sound_s = sound_s / max(sound_s);
        sound_s = sound_s * x_range + min(x_avg); 
    end
    h_sound = plot(ax, t_s, sound_s, 'linew', linew); 
    h_sound.Color = col_sound; 
end

% plot windows of interest 
if ~isempty(win_start_times) && ~isempty(win_end_times)
    plot_erp_windows(ax, x_avg, fs,...
        win_start_times, win_end_times, col_win);
end

% if multiple epochs or subjects were passed, get mean and CIs
if size(x, 1) > 1 && plot_ci
    if strcmp(ci_method, 'boot')
        % bootstrap confidence intervals [0.7,0.7,0.7]
        [ci_high, ci_low] = bootstrapCI(x, par.nBoot); 
    elseif strcmp(ci_method, 'analytical')
        ci = std(x, [], 1) / sqrt(size(x, 1)) * norminv(1-0.05);
        ci_high = x_avg + ci; 
        ci_low = x_avg - ci;
    else 
        error('CI method "%s" not implemented', ci_method);
    end
    t_ci = [t, fliplr(t)]; 
    ci = [ci_high, fliplr(ci_low)]; 
    % plot CIs as shaded regions
    fill(ax, t_ci, ci, col, ...
        'facealpha',0.3, ...
        'LineStyle','none')
else
    ci = mean(x,1); 
end

h = plot(ax, t, x_avg,...
    'color',col,...
    'linewidth',linew); 

ax.XAxis.Visible = 'off'; 
ax.YAxis.Visible = 'on'; 
ax.TickDir = 'out'; 
ax.XLim = [t(1), t(end)]; 

if isnumeric(y_max) 
    % mirror y-axs limits if needed 
    if length(y_max)==1
        y_max = [-y_max,y_max]; 
    end
    ax.YLim = [min(min(ci),y_max(1)), max(max(ci),y_max(2))]; 
    if ~any(isinf(y_max))
        ax.YTick = y_max; 
    end
else
    ax.YTick = [min(ax.YTick), max(ax.YTick)];
end

ax.FontSize = fontsize; 

