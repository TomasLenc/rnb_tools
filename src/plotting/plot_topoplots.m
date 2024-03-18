function [f, ax, cbar] = plot_topoplots(vals, chanlocs, varargin)
% Plot an array of topopolots of any set of values estimated at each
% electrode, separately for each condition. 
% 
% Parameters
% ----------
% vals : array_like, shape=[n_electrodes, n_conditions]
%     Two dimensional array with values to be plotted across electrodes on
%     one dimension, and conditions on the other dimension. One topoplot
%     will be plotted per condition into separate axes in teh figure. 
% chanlocs : struct
%     eeglab style chanlocs, that will be passed to the topoplot eeglab
%     function. 
% lab : string, optional
%     Name of the plotted variable, will be shown as the label on the
%     colorbar. 
% mark_chan_idx : array of integers, optional
%     Indices of channels that will be marked with a point in the figure. 
% cmap : array_like, shape=[n_colors, 3]
%     Array 
% mark_chan_idx : array of integers, optional
%     Indices of channels that will be marked with a point in the figure. 
% 
% Returns 
% -------
% f : figure object 
%     Handle to the generated figure. 
% ax : axes object 
%     Handle to the generated axes. 
% cbar : colobar object 
%     Handle to the generated colobar. 
% 

parser = inputParser(); 

addParameter(parser, 'lab', ''); 
addParameter(parser, 'mark_chan_idx', []); 
addParameter(parser, 'cmap', parula); 
addParameter(parser, 'ax', []); 

parse(parser, varargin{:}); 

lab = parser.Results.lab; 
mark_chan_idx = parser.Results.mark_chan_idx; 
cmap = parser.Results.cmap; 
ax = parser.Results.ax; 


if isempty(mark_chan_idx)
    show_elec = 'off'; 
else
    show_elec = 'on'; 
end

cmax = 1.0 * prctile(vals(:), 100); 

% find which dimension 
idx_dim_chans = find(size(vals) == length(chanlocs)); 
dims = [1 : ndims(vals)]; 
n_cond = size(vals, dims(dims ~= idx_dim_chans)); 

if isempty(ax)    
    
    if n_cond == 1
        fig_size = [565 604 180 146] / 1.3; 
        marker_size = 4; 
    else 
        fig_size = [565 604 100*n_cond 146] / 1.3; 
        marker_size = 2; 
    end
    
    f = figure('color', 'white', 'Position', fig_size); 
    pnl = panel(f); 
    pnl.pack('h', n_cond);     
    ax = nan(1, n_cond); 
    for i_cond=1:n_cond
        ax(i_cond) = pnl(i_cond).select(); 
    end
else 
    f = []; 
end

for i_cond=1:n_cond

    axes(ax(i_cond)); 
    
    idx = {i_cond, i_cond}; 
    idx{idx_dim_chans} = ':'; 
    
    topoplot(vals(idx{:}), chanlocs, ...
             'style', 'map', ...
             'maplimits', [0, cmax], ...
             'gridscale', 256, ...
             'colormap', cmap, ...
             'electrodes', show_elec, ...
             'emarker2', {mark_chan_idx,'o','r',marker_size});

    if i_cond == n_cond
        cbar = colorbar(); 
        cbar.Ticks = [0, cmax]; 
        cbar.TickLabels = {0, floor(cmax*100)/100}; 
        cbar.FontSize = 12; 
        cbar.Label.String = lab; 
        cbar.Label.Position(1) = cbar.Label.Position(1) - cbar.Label.Position(1)*0.5; 
    end
    
    if n_cond > 1
        title(ax(i_cond), sprintf('%d', i_cond)); 
    end
    
end

pnl.de.margin = [1, 1, 1, 1]; 
pnl.margin = [1, 1, 20, 1]; 
