function ax = plot_circ(ph, varargin)
% This function plots timeseries defined as a vector of phase angles on a
% circular plot. Optionally also plots the mean vector. 
% 
% Parameters
% ----------
% ph : array_like
%     Vector with phase angle assigned to each onset time. 
% r : float, optional 
%     Mean vector length. 
% mean_ph : float, optional 
%     Phase angle of the mean vector. 
% ax : axes, optional 
%     Instance of axes to plot into. 
% mean_ph : float, optional 
%     Phase angle of the mean vector. 
% linew : float, optional 
%     Width of the line plotting mean vector. 
% fontsize : float, optional 
%     Fontsize.
% col_ind : array of [int int int], optional 
%     RGB color for the individual points. 
% col_mean : array of [int int int], optional 
%     RGB color for the mean vector.
% 
% Returns 
% -------
% closest_period : axes
%   Instance of axes. 
% 

parser = inputParser; 

addParameter(parser, 'ax', []); 
addParameter(parser, 'mean_ph', []); 
addParameter(parser, 'r', []); 
addParameter(parser, 'linew', 1.7); 
addParameter(parser, 'fontsize', 12); 
addParameter(parser, 'col_ind', [115, 115, 115]/255); 
addParameter(parser, 'col_mean', [196, 70, 16]/255); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
mean_ph = parser.Results.mean_ph; 
r = parser.Results.r; 
linew = parser.Results.linew; 
fontsize = parser.Results.fontsize; 
col_ind = parser.Results.col_ind; 
col_mean = parser.Results.col_mean; 

if isempty(ax)
    f = figure('color', 'white', 'position', [680 838 153 124]);
    ax = polaraxes;
end

polarplot(ax, ph, ones(size(ph)), 'o', 'color', col_ind);
hold(ax, 'on');
if ~isempty(mean_ph) && ~isempty(r)
    polarplot(ax, [mean_ph, mean_ph]', [0, r]', ...
              'color', col_mean, 'linew', linew);
end

ax.ThetaAxisUnits = 'radians';
ax.ThetaTick = [0, pi/2, pi, 3*pi/2];
ax.ThetaTickLabel = [];
ax.RTick = [];
ax.RLim = [0, 1.1];
ax.FontSize = fontsize;

